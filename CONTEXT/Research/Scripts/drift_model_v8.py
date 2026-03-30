#!/usr/bin/env python3
"""
Дрейф ЭВО v8: амплитуда A из нелинейного КУБК
===============================================
B (Математика), 2026-03-29 | A-DIR-AMP-002
 
Амплитуда: из условия ⟨p_канал(ψ)·|cos_tilt(ψ)|⟩_ψ = p_СЭ.
  - p_канал(ψ) = p₀/(1−a·sinψ)² (компрессия канала)
  - cos_tilt(ψ) = проекция деформированной нормали на радиальную
  - Всё — через статистику столкновений, без макро-сущностей.
 
v6 (channel_forces, chaos_forces, compute_drift) — без изменений.
Новое: find_A_eq (модуль определения амплитуды).
"""
 
import numpy as np
from scipy.optimize import brentq
 
d_A = 1.0; v_0 = 1.0; m_a = 1.0
n_a = 16/(27*np.sqrt(2)*d_A**3)
p_SE = m_a*v_0**2/(4*d_A**2)
G_GAS = 1.0/3.0; Q_MIN = 3
R_major = 2500.0; r_minor = 50.0
N_TH = 512; N_PSI = 1024  # высокое разрешение по ψ для точного интеграла
 
 
def wave_params(r, m_S=1):
    k=1.0/(m_S*r); kz=k/np.sqrt(2)
    return {'k_z':kz,'Lz_min':np.sqrt(2)*np.pi*m_S*r}
 
def close_params(v_n_over_v0=5.0, v_ratio=2.0):
    v_n=v_n_over_v0*v_0; v_ch=v_n/v_ratio
    n_ch=v_0**2/(4*d_A**2*G_GAS*v_ch**2)
    return {'v_n':v_n,'v_ch':v_ch,'n_ch':n_ch,'n_ratio':n_ch/n_a}
 
 
# ============================================================
# НЕЛИНЕЙНЫЙ КУБК: ОПРЕДЕЛЕНИЕ АМПЛИТУДЫ
# ============================================================
def kubk_residual(A, r, n_ch0, v_ch, m_S=1):
    """
    Невязка КУБК: ⟨p_канал(ψ)·|cos_tilt(ψ)|⟩_ψ − p_СЭ.
 
    p_канал(ψ) = G_gas · n_ch(ψ) · m_a · v_ch²
    n_ch(ψ) = n_ch0 · r² / (r − A·sinψ)²  (компрессия)
 
    cos_tilt(ψ): проекция деформированной нормали на радиальное
    направление n̂₀. Для волны δρ = −A·sin(ψ):
      Наклон: dρ/ds = −A·k_z·cos(ψ) вдоль z, −A·(k_θ/r)·cos(ψ) вдоль θ.
      |cos_tilt| = 1/√(1 + (A·k_z·cos(ψ))²·2) для 45° режима.
      При малых Ak: |cos_tilt| ≈ 1 − (Ak)²cos²ψ.
    """
    if A <= 0: return -p_SE  # при A=0: p_ch·1 = p_ch0 = p_SE → невязка=0
    if A >= r: return np.inf
 
    a = A / r
    wp = wave_params(r, m_S)
    kz = wp['k_z']
    ak = A * kz
 
    psi = np.linspace(0, 2*np.pi, N_PSI, endpoint=False)
    dpsi = 2*np.pi / N_PSI
 
    # Компрессия канала
    denom = 1 - a * np.sin(psi)
    # Защита от деления на ноль
    denom = np.maximum(denom, 0.01)
    n_ch_local = n_ch0 / denom**2
 
    # Давление канала (газокинетическое)
    p_local = G_GAS * n_ch_local * m_a * v_ch**2
 
    # Наклон нормали: для волны δρ = −A·sin(k_z·z + k_θ·θ − ωt)
    # Градиент фазы в z: k_z, в θ: k_θ/r = k_z (для 45°)
    # Полный наклон: |∇δρ| = A·|k⃗|·|cosψ| = A·√2·k_z·|cosψ|
    # cos_tilt = 1/√(1 + (∇δρ)²) = 1/√(1 + 2·(Ak_z)²·cos²ψ)
    grad_sq = 2 * (ak * np.cos(psi))**2
    cos_tilt = 1.0 / np.sqrt(1 + grad_sq)
 
    # Среднее эффективное давление
    p_eff_mean = np.sum(p_local * cos_tilt) * dpsi / (2*np.pi)
 
    return p_eff_mean - p_SE
 
 
def find_A_eq(r, n_ch0, v_ch, m_S=1):
    """
    Найти A_eq из ⟨p_канал(ψ)·|cos_tilt(ψ)|⟩ = p_СЭ.
    """
    # Проверяем знак невязки при A→0+ и A→r-
    res_small = kubk_residual(0.01*d_A, r, n_ch0, v_ch, m_S)
    res_big = kubk_residual(0.95*r, r, n_ch0, v_ch, m_S)
 
    if res_small >= 0 and res_big >= 0:
        # Невязка всегда ≥ 0 → p_ch всегда ≥ p_SE → нет решения кроме A=0
        # Но проверим: может быть невязка при A=0 ровно 0 (базовый КУБК)
        return 0.0, False, "p_ch ≥ p_SE для всех A"
 
    if res_small * res_big > 0:
        # Одного знака → ищем минимум невязки
        from scipy.optimize import minimize_scalar
        result = minimize_scalar(lambda A: abs(kubk_residual(A, r, n_ch0, v_ch, m_S)),
                                 bounds=(0.01, 0.95*r), method='bounded')
        if abs(kubk_residual(result.x, r, n_ch0, v_ch, m_S)) < 1e-6 * p_SE:
            return result.x, True, "минимум невязки"
        return 0.0, False, f"нет корня (res_small={res_small:.3e}, res_big={res_big:.3e})"
 
    # Разные знаки → есть корень
    try:
        A_eq = brentq(lambda A: kubk_residual(A, r, n_ch0, v_ch, m_S),
                      0.01*d_A, 0.95*r, xtol=1e-6)
        return A_eq, True, "brentq"
    except:
        return 0.0, False, "brentq failed"
 
 
# ============================================================
# v6 КОМПОНЕНТЫ (без изменений)
# ============================================================
def deformed_normal(TH,PSI,chi,s_sh,kz,A):
    c=A*s_sh*kz*np.cos(PSI)
    n0x=np.cos(TH);n0y=np.sin(TH);ex=-np.sin(TH);ey=np.cos(TH)
    nx=n0x-chi*c*ex;ny=n0y-chi*c*ey;nz=-c*np.ones_like(TH)
    mag=np.sqrt(nx**2+ny**2+nz**2); return nx/mag,ny/mag,nz/mag
 
def channel_forces(TH,PSI,chi,s_sh,v_n,v_ch,n_ch,kz,A):
    nx,ny,nz=deformed_normal(TH,PSI,chi,s_sh,kz,A)
    ex=-np.sin(TH);ey=np.cos(TH)
    vnx=v_n/np.sqrt(2)*chi*ex;vny=v_n/np.sqrt(2)*chi*ey
    vnz=v_n/np.sqrt(2)*np.ones_like(TH)
    cos_tilt=(vnx*nx+vny*ny+vnz*nz)/v_n
    p_base=G_GAS*n_ch*m_a*v_ch**2
    Fz_b=p_base*(-nz);Fth_b=p_base*(-(-np.sin(TH)*nx+np.cos(TH)*ny))
    drift=G_GAS*n_ch*m_a*v_ch**2*cos_tilt**2*(cos_tilt<0).astype(float)
    return Fz_b+drift*nz,Fth_b+drift*(-np.sin(TH)*nx+np.cos(TH)*ny),p_base
 
def chaos_forces(TH,PSI,chi,s_sh,kz,A):
    nx,ny,nz=deformed_normal(TH,PSI,chi,s_sh,kz,A)
    return -p_SE*nz,-p_SE*(-np.sin(TH)*nx+np.cos(TH)*ny)
 
def compute_drift(chi,s_sh,v_n,v_ch,n_ch,A,m_S=1,q=3.0):
    assert q>=Q_MIN
    wp=wave_params(r_minor,m_S);kz=wp['k_z']
    th=np.linspace(0,2*np.pi,N_TH,endpoint=False)
    psi=np.linspace(0,2*np.pi,N_PSI,endpoint=False)
    TH,PSI=np.meshgrid(th,psi,indexing='ij')
    dth=2*np.pi/N_TH;dpsi=2*np.pi/N_PSI
    avg=lambda F:np.sum(F)*dth*dpsi/(2*np.pi)
    Fz_ch,Fth_ch,p_ch=channel_forces(TH,PSI,chi,s_sh,v_n,v_ch,n_ch,kz,A)
    Fz_se,Fth_se=chaos_forces(TH,PSI,chi,s_sh,kz,A)
    Fth_avg=avg(Fth_ch)+avg(Fth_se)
    gamma=n_a*m_a*v_0*(2*np.pi*r_minor)
    Omega=-Fth_avg/(gamma*r_minor) if gamma*r_minor!=0 else 0.0
    coupling=float(s_sh)*float(chi)*q/np.sqrt(1+q**2) if s_sh!=0 else 0.0
    return coupling*Omega*r_minor, {'Omega':Omega,'Fth_ch':avg(Fth_ch),
        'Fth_se':avg(Fth_se),'p_ch':float(p_ch),'p_ratio':float(p_ch)/p_SE}
 
 
# ============================================================
# MAIN
# ============================================================
def main():
    H="="*90
    print("╔═══════════════════════════════════════════════════════════════════════╗")
    print("║  DRIFT v8: АМПЛИТУДА A ИЗ НЕЛИНЕЙНОГО КУБК                        ║")
    print("║  A-DIR-AMP-002 | ⟨p_канал·|cos_tilt|⟩ = p_СЭ                       ║")
    print("╚═══════════════════════════════════════════════════════════════════════╝")
 
    # --- Диагностика невязки ---
    print(f"\n{H}")
    print("1. ДИАГНОСТИКА НЕВЯЗКИ КУБК(A)")
    print(H)
 
    ch = close_params(v_n_over_v0=5.0)
    v_n=ch['v_n'];v_ch=ch['v_ch'];n_ch=ch['n_ch']
 
    print(f"\n  |v_n|={v_n:.1f}v₀, |v|={v_ch:.2f}v₀, n_ch/n_a={ch['n_ratio']:.3f}")
    print(f"  p_канал₀ = {G_GAS*n_ch*m_a*v_ch**2:.4f}, p_СЭ = {p_SE:.4f}")
 
    print(f"\n  {'A [d_A]':>10} {'a=A/r':>8} {'⟨p·|cos|⟩':>12} {'p_СЭ':>8} {'невязка':>12}")
    print(f"  {'─'*55}")
    for A_val in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 7.07, 10.0, 15.0, 20.0, 30.0, 40.0]:
        if A_val >= r_minor: continue
        res = kubk_residual(A_val, r_minor, n_ch, v_ch)
        p_eff = res + p_SE
        print(f"  {A_val:10.2f} {A_val/r_minor:8.4f} {p_eff:12.6f} {p_SE:8.4f} {res:+12.6e}")
 
    # --- Поиск A_eq ---
    print(f"\n{H}")
    print("2. ПОИСК A_eq")
    print(H)
 
    A_eq, found, method = find_A_eq(r_minor, n_ch, v_ch)
    if found:
        a_eq = A_eq/r_minor
        kz = wave_params(r_minor)['k_z']
        print(f"\n  A_eq = {A_eq:.4f} d_A")
        print(f"  A_eq/r = {a_eq:.6f}")
        print(f"  A_eq·k_z = {A_eq*kz:.6f}")
        print(f"  Метод: {method}")
        print(f"  Проверка: невязка = {kubk_residual(A_eq, r_minor, n_ch, v_ch):.2e}")
 
        # Проверка n_ch(ψ) < n_a
        a = A_eq/r_minor
        n_max = n_ch / (1-a)**2 if a < 1 else np.inf
        print(f"  n_ch_max (на вершине волны) = {n_max:.4f}/d³ = {n_max/n_a:.3f}·n_a "
              f"{'✅' if n_max<n_a else '⚠️ > n_a'}")
    else:
        print(f"\n  A_eq не найден: {method}")
        print(f"  Используем A = 0 (гладкая ГП) или A как параметр.")
        A_eq = 7.07  # fallback
 
    # --- Сканы ---
    print(f"\n{H}")
    print("3. СКАНЫ A_eq ПО ПАРАМЕТРАМ")
    print(H)
 
    print(f"\n  3a) A_eq(|v_n|):")
    print(f"  {'|v_n|/v₀':>10} {'n_ch/n_a':>9} {'A_eq':>8} {'a=A/r':>8} {'A·k_z':>8} {'found':>6}")
    for vr in [3.0,4.0,5.0,7.0,10.0]:
        c2=close_params(v_n_over_v0=vr)
        A2,ok,_=find_A_eq(r_minor,c2['n_ch'],c2['v_ch'])
        kz=wave_params(r_minor)['k_z']
        if ok:
            print(f"  {vr:10.1f} {c2['n_ratio']:9.3f} {A2:8.4f} {A2/r_minor:8.6f} {A2*kz:8.6f} {'✅':>6}")
        else:
            print(f"  {vr:10.1f} {c2['n_ratio']:9.3f} {'—':>8} {'—':>8} {'—':>8} {'❌':>6}")
 
    print(f"\n  3b) A_eq(r):")
    for r_val in [25.0,50.0,100.0]:
        c2=close_params(v_n_over_v0=5.0)
        A2,ok,_=find_A_eq(r_val,c2['n_ch'],c2['v_ch'])
        kz2=wave_params(r_val)['k_z']
        if ok:
            print(f"  r={r_val:5.0f}: A_eq={A2:.4f}, a={A2/r_val:.6f}, Ak={A2*kz2:.6f}")
 
    # --- V_∥ при A_eq ---
    if found and A_eq > 0:
        print(f"\n{H}")
        print("4. V_∥ ПРИ САМОСОГЛАСОВАННОМ A")
        print(H)
 
        q = 3.0
        print(f"\n  {'χ':>3} {'s':>3} {'V_∥':>13} {'Ω':>13}")
        print(f"  {'─'*40}")
        R = {}
        for chi in [+1,-1]:
            for s in [+1,0,-1]:
                V,d = compute_drift(chi,s,v_n,v_ch,n_ch,A_eq,q=q)
                R[(chi,s)] = (V,d)
                print(f"  {chi:+3d} {s:+3d} {V:+13.4e} {d['Omega']:+13.4e}")
 
        # Симметрии
        print(f"\n  Симметрии:")
        for s in [+1,-1]:
            Vp=R[(+1,s)][0];Vm=R[(-1,s)][0]
            dv=abs(Vp-Vm)/max(abs(Vp),1e-30)
            print(f"    V_∥ χ-чётно (s={s:+d}): δ={dv:.2e} → {'✅' if dv<1e-6 else '❌'}")
        for chi in [+1,-1]:
            Vf=R[(chi,+1)][0];Vr=R[(chi,-1)][0]
            dr=abs(Vf+Vr)/max(abs(Vf),1e-30)
            print(f"    Реверс (χ={chi:+d}): δ={dr:.2e} → {'✅' if dr<1e-6 else '❌'}")
 
        V_ref = abs(R[(+1,+1)][0])
        print(f"\n  V_∥ = {V_ref:.6e} · v₀")
        print(f"  V_∥/v_s = {V_ref/(3*np.sqrt(2)*v_0):.6e}")
 
        # Сравнение
        V_old,_ = compute_drift(+1,+1,v_n,v_ch,n_ch,7.07,q=q)
        print(f"\n  Сравнение: V(A=7.07) = {abs(V_old):.4e}, V(A_eq={A_eq:.2f}) = {V_ref:.4e}")
 
    print(f"\n{H}")
 
 
if __name__=='__main__':
    main()
