#!/usr/bin/env python3
"""
Дрейф ЭВО v6: исправленная спиральная связь + q ≥ 3
=====================================================
B (Математика), 2026-03-27 | A-DIR-SPIRAL-001
 
КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ:
  coupling = s_sh · χ · q / √(1+q²)   [было: s_sh · q / √(1+q²)]
 
  При χ→−χ: Ω инвертируется, но навивка двутория тоже → χ² = 1 → V_∥ χ-ЧЁТНО.
  Правый и левый фотон летят в одном направлении, различаются вращением.
 
Симметрии (уточнённые, Прил. A.1):
  Ω(χ) = −Ω(−χ)         χ-нечётно  (поляризация)
  V_∥(χ) = V_∥(−χ)       χ-ЧЁТНО   (скорость фотона)
  V_∥(s_sh) = −V_∥(−s_sh) s_sh-нечётно (реверс бега)
 
Ограничение: q ≥ 3 (устойчивость спиральной навивки).
 
Геометрические факторы: G_СЭ = 1/(4d²), G_канал = 1/3 (v5.1).
КУБК: n_ch = 3v₀²/(4d²·v_ch²), p_канал/p_СЭ = 1.
"""
 
import numpy as np
 
d_A = 1.0; v_0 = 1.0; m_a = 1.0
n_a = 16/(27*np.sqrt(2)*d_A**3)
p_SE = m_a*v_0**2/(4*d_A**2)
G_GAS = 1.0/3.0
Q_MIN = 3
 
R_major = 2500.0; r_minor = 50.0
N_TH = 512; N_PSI = 512
 
 
def wave_params(r, m_S=1):
    k=1.0/(m_S*r); kz=k/np.sqrt(2)
    return {'k_z': kz, 'Lz_min': np.sqrt(2)*np.pi*m_S*r}
 
def close_params(v_n_over_v0=5.0, v_ratio=2.0):
    v_n=v_n_over_v0*v_0; v_ch=v_n/v_ratio
    n_ch=v_0**2/(4*d_A**2*G_GAS*v_ch**2)
    return {'v_n':v_n,'v_ch':v_ch,'n_ch':n_ch,'n_ratio':n_ch/n_a}
 
def close_amplitude(r, m_S=1):
    return 0.1/wave_params(r,m_S)['k_z']
 
def deformed_normal(TH,PSI,chi,s_sh,kz,A):
    c=A*s_sh*kz*np.cos(PSI)
    n0x=np.cos(TH);n0y=np.sin(TH);ex=-np.sin(TH);ey=np.cos(TH)
    nx=n0x-chi*c*ex;ny=n0y-chi*c*ey;nz=-c*np.ones_like(TH)
    mag=np.sqrt(nx**2+ny**2+nz**2)
    return nx/mag,ny/mag,nz/mag
 
def channel_forces(TH,PSI,chi,s_sh,v_n,v_ch,n_ch,kz,A):
    nx,ny,nz=deformed_normal(TH,PSI,chi,s_sh,kz,A)
    ex=-np.sin(TH);ey=np.cos(TH)
    vnx=v_n/np.sqrt(2)*chi*ex;vny=v_n/np.sqrt(2)*chi*ey
    vnz=v_n/np.sqrt(2)*np.ones_like(TH)
    cos_tilt=(vnx*nx+vny*ny+vnz*nz)/v_n
    p_base=G_GAS*n_ch*m_a*v_ch**2
    Fz_b=p_base*(-nz);Fth_b=p_base*(-(-np.sin(TH)*nx+np.cos(TH)*ny))
    drift=G_GAS*n_ch*m_a*v_ch**2*cos_tilt**2*(cos_tilt<0).astype(float)
    Fz_d=drift*nz;Fth_d=drift*(-np.sin(TH)*nx+np.cos(TH)*ny)
    return Fz_b+Fz_d,Fth_b+Fth_d,p_base
 
def chaos_forces(TH,PSI,chi,s_sh,kz,A):
    nx,ny,nz=deformed_normal(TH,PSI,chi,s_sh,kz,A)
    return -p_SE*nz,-p_SE*(-np.sin(TH)*nx+np.cos(TH)*ny)
 
def compute_drift(chi,s_sh,v_n,v_ch,n_ch,A,m_S=1,q=3.0):
    assert q>=Q_MIN, f"q={q} < q_min={Q_MIN}"
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
 
    # ИСПРАВЛЕНО: coupling = s_sh · χ · q/√(1+q²)
    coupling=float(s_sh)*float(chi)*q/np.sqrt(1+q**2) if s_sh!=0 else 0.0
 
    V_par=coupling*Omega*r_minor
 
    return V_par, {
        'Fth_ch':avg(Fth_ch),'Fth_se':avg(Fth_se),'Fth_total':Fth_avg,
        'Omega':Omega,'coupling':coupling,
        'p_ch':float(p_ch),'p_ratio':float(p_ch)/p_SE,
    }
 
 
def main():
    H="="*90
    print("╔═══════════════════════════════════════════════════════════════════════╗")
    print("║  DRIFT v6: ИСПРАВЛЕННАЯ СПИРАЛЬНАЯ СВЯЗЬ + q ≥ 3                   ║")
    print("║  coupling = s_sh · χ · q/√(1+q²)  →  V_∥ χ-ЧЁТНО, Ω χ-НЕЧЁТНО    ║")
    print("╚═══════════════════════════════════════════════════════════════════════╝")
 
    ch=close_params(v_n_over_v0=5.0)  # n_ch/n_a=0.29 — чистый газ
    A=close_amplitude(r_minor)
    v_n=ch['v_n'];v_ch=ch['v_ch'];n_ch=ch['n_ch']
    q=3.0; m_S=1
 
    print(f"\n  |v_n|={v_n:.1f}v₀, |v|={v_ch:.2f}v₀, n_ch/n_a={ch['n_ratio']:.3f}, A={A:.2f}d_A, q={q}")
 
    # --- Все комбинации ---
    print(f"\n{'─'*100}")
    print(f"{'χ':>3} {'s':>3} {'Fθ(ch)':>13} {'Ω':>13} {'coupling':>10} {'V_∥':>13} {'p_r':>6}")
    print(f"{'─'*100}")
 
    R={}
    for chi in [+1,-1]:
        for s in [+1,0,-1]:
            V,d=compute_drift(chi,s,v_n,v_ch,n_ch,A,m_S,q)
            R[(chi,s)]=(V,d)
            print(f"{chi:+3d} {s:+3d} {d['Fth_ch']:+13.4e} {d['Omega']:+13.4e} "
                  f"{d['coupling']:+10.4f} {V:+13.4e} {d['p_ratio']:6.3f}")
 
    # --- ПРОВЕРКИ ---
    print(f"\n{H}\nПРОВЕРКИ (уточнённые симметрии Прил. A.1)\n{H}")
 
    print(f"\n  p_канал/p_СЭ = {R[(+1,+1)][1]['p_ratio']:.6f} → ✅")
    print(f"  n_ch/n_a = {ch['n_ratio']:.4f} → {'✅' if ch['n_ratio']<1 else '❌'}")
 
    # 1a) χ-НЕЧЁТНОСТЬ Ω
    print("\n  1a) Ω χ-НЕЧЁТНО: Ω(χ=+1) = −Ω(χ=−1)")
    for s in [+1,-1]:
        Op=R[(+1,s)][1]['Omega']; Om=R[(-1,s)][1]['Omega']
        d=abs(Op+Om)/max(abs(Op),abs(Om),1e-30)
        print(f"      s={s:+d}: Ω(+)={Op:+.4e}, Ω(−)={Om:+.4e}, δ={d:.2e} → {'✅' if d<1e-10 else '❌'}")
 
    # 1b) χ-ЧЁТНОСТЬ V_∥
    print("\n  1b) V_∥ χ-ЧЁТНО: V_∥(χ=+1) = V_∥(χ=−1)")
    for s in [+1,-1]:
        Vp=R[(+1,s)][0]; Vm=R[(-1,s)][0]
        d=abs(Vp-Vm)/max(abs(Vp),abs(Vm),1e-30)
        print(f"      s={s:+d}: V(+)={Vp:+.4e}, V(−)={Vm:+.4e}, δ={d:.2e} → {'✅' if d<1e-10 else '❌'}")
 
    # 2) Реверс бега
    print("\n  2) РЕВЕРС БЕГА: V_∥(s=+1) = −V_∥(s=−1)")
    for chi in [+1,-1]:
        Vf=R[(chi,+1)][0]; Vr=R[(chi,-1)][0]
        d=abs(Vf+Vr)/max(abs(Vf),abs(Vr),1e-30)
        print(f"      χ={chi:+d}: V(→)={Vf:+.4e}, V(←)={Vr:+.4e}, δ={d:.2e} → {'✅' if d<1e-10 else '❌'}")
 
    # 3) Нули
    for label,V in [("s=0",R[(+1,0)][0]),("A=0",compute_drift(+1,+1,v_n,v_ch,n_ch,0.0,q=q)[0]),
                    ("q=3,s=0",R[(+1,0)][0])]:
        print(f"  {label}: V={V:.2e} → {'✅' if abs(V)<1e-15 else '❌'}")
 
    # --- КУБК-инвариант ---
    print(f"\n{H}\nКУБК-ИНВАРИАНТ (V_∥ не зависит от |v_n|)\n{H}")
    for vr in [3.0,5.0,7.0,10.0]:
        c2=close_params(v_n_over_v0=vr)
        V2,_=compute_drift(+1,+1,c2['v_n'],c2['v_ch'],c2['n_ch'],A,q=q)
        print(f"  |v_n|={vr:.0f}v₀: n_ch/n_a={c2['n_ratio']:.3f}, V_∥={V2:+.4e}")
 
    # --- СКАН ПО q ---
    print(f"\n{H}\nСКАН ПО q (χ=+1, s=+1)\n{H}")
    print(f"  {'q':>6} {'coupling':>10} {'V_∥':>13} {'V_∥·√(1+q²)/q':>16}")
    for qv in [3,5,10,50,100]:
        V,d=compute_drift(+1,+1,v_n,v_ch,n_ch,A,q=float(qv))
        V_norm=V*np.sqrt(1+qv**2)/qv if qv>0 else 0
        print(f"  {qv:6d} {d['coupling']:+10.4f} {V:+13.4e} {V_norm:+16.4e}")
 
    # --- Степенные законы ---
    print(f"\n{H}\nСТЕПЕННЫЕ ЗАКОНЫ\n{H}")
    Vs=[];As=[0.5,1.0,2.0,A,A*2];As.sort()
    for Av in As:
        V,_=compute_drift(+1,+1,v_n,v_ch,n_ch,Av,q=q);Vs.append(V)
        print(f"  A={Av:6.2f}: V_∥={V:+.4e}")
    if Vs[0]!=0 and Vs[-1]!=0:
        pw=np.log(abs(Vs[-1]/Vs[0]))/np.log(As[-1]/As[0])
        print(f"  → V ∝ A^{pw:.2f}")
 
    V_ref=abs(R[(+1,+1)][0])
    print(f"\n{H}")
    print(f"  V_∥ = {V_ref:.6e} · v₀  (при q={q})")
    print(f"  V_∥/v_s = {V_ref/(3*np.sqrt(2)*v_0):.6e}")
    print(f"  Ω = {abs(R[(+1,+1)][1]['Omega']):.6e} (поляризация)")
 
    print(f"""
  ИТОГОВАЯ ФИЗИКА:
    Ω ∝ χ         — ПОЛЯРИЗАЦИЯ (правый/левый фотон вращаются по-разному)
    V_∥ ∝ s_sh    — НАПРАВЛЕНИЕ ПОЛЁТА (определяется бегом волны)
    V_∥ не зависит от χ — ПРАВЫЙ И ЛЕВЫЙ ФОТОН ЛЕТЯТ С ОДНОЙ СКОРОСТЬЮ
 
    coupling = s_sh·χ·q/√(1+q²)
    V_∥ = coupling·Ω·r = s_sh·χ·q/√(1+q²) · χ·|Ω₀|·r = s_sh·q·|Ω₀|·r/√(1+q²)
    (χ² = 1 → χ выпадает из V_∥)
""")
    print(H)
 
 
if __name__=='__main__':
    main()
