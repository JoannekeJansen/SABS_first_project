
"""This module contains a Paci2013_atrial cardiac cell model

The module was autogenerated from a gotran ode file
"""
from __future__ import division
from collections import OrderedDict
import ufl

from cbcbeat.dolfinimport import *
from cbcbeat.cellmodels import CardiacCellModel

class Paci2013_atrial(CardiacCellModel):
    def __init__(self, params=None, init_conditions=None):
        """
        Create cardiac cell model

        *Arguments*
         params (dict, :py:class:`dolfin.Mesh`, optional)
           optional model parameters
         init_conditions (dict, :py:class:`dolfin.Mesh`, optional)
           optional initial conditions
        """
        CardiacCellModel.__init__(self, params, init_conditions)

    @staticmethod
    def default_parameters():
        "Set-up and return default parameters."
        params = OrderedDict([("PkNa", 0.03),
                              ("g_Na", 6646.185),
                              ("g_CaL", 8.635702e-05),
                              ("tau_fCa", 0.002),
                              ("g_Kr", 29.8667),
                              ("L0", 0.025),
                              ("Q", 2.3),
                              ("g_Ks", 2.041),
                              ("g_K1", 19.1925),
                              ("E_f", -0.017),
                              ("g_f", 30.10312),
                              ("g_b_Na", 0.9),
                              ("g_b_Ca", 0.69264),
                              ("Km_K", 1.0),
                              ("Km_Na", 40.0),
                              ("PNaK", 1.4731392),
                              ("KmCa", 1.38),
                              ("KmNai", 87.5),
                              ("Ksat", 0.1),
                              ("alpha", 2.8571432),
                              ("gamma", 0.35),
                              ("kNaCa", 2450.0),
                              ("KPCa", 0.0005),
                              ("g_PCa", 0.4125),
                              ("g_to", 59.8077),
                              ("Cao", 1.8),
                              ("Cm", 7.86671e-11),
                              ("F", 96485.3415),
                              ("Ki", 150.0),
                              ("Ko", 5.4),
                              ("Nao", 151.0),
                              ("R", 8.314472),
                              ("T", 310.0),
                              ("V_SR", 465.2),
                              ("Vc", 7012.0),
                              ("Buf_C", 0.25),
                              ("Buf_SR", 10.0),
                              ("Kbuf_C", 0.001),
                              ("Kbuf_SR", 0.3),
                              ("Kup", 0.00025),
                              ("V_leak", 0.00044444),
                              ("VmaxUp", 0.22),
                              ("a_rel", 16.464),
                              ("b_rel", 0.25),
                              ("c_rel", 8.232),
                              ("tau_g", 0.002)])
        return params

    @staticmethod
    def default_initial_conditions():
        "Set-up and return default initial conditions."
        ic = OrderedDict([("V", -0.0687338234522),
                          ("m", 0.141183142078),
                          ("h", 0.642108593995),
                          ("j", 0.173566329483),
                          ("d", 0.000127632520742),
                          ("f1", 0.980384004336),
                          ("f2", 0.99995300671),
                          ("fCa", 0.997346890769),
                          ("Xr1", 0.0257889110986),
                          ("Xr2", 0.40504667874),
                          ("Xs", 0.0447460799149),
                          ("Xf", 0.0607988713875),
                          ("q", 0.776163826643),
                          ("r", 0.000503296941001),
                          ("Nai", 14.4424010544),
                          ("g", 1.0),
                          ("Cai", 4.49232909235e-05),
                          ("Ca_SR", 0.149980051222)])
        return ic

    def _I(self, v, s, time):
        """
        Original gotran transmembrane current dV/dt
        """
        time = time if time else Constant(0.0)

        # Assign states
        V = v
        assert(len(s) == 17)
        m, h, j, d, f1, f2, fCa, Xr1, Xr2, Xs, Xf, q, r, Nai, g, Cai, Ca_SR = s

        # Assign parameters
        PkNa = self._parameters["PkNa"]
        g_Na = self._parameters["g_Na"]
        g_CaL = self._parameters["g_CaL"]
        g_Kr = self._parameters["g_Kr"]
        g_Ks = self._parameters["g_Ks"]
        g_K1 = self._parameters["g_K1"]
        E_f = self._parameters["E_f"]
        g_f = self._parameters["g_f"]
        g_b_Na = self._parameters["g_b_Na"]
        g_b_Ca = self._parameters["g_b_Ca"]
        Km_K = self._parameters["Km_K"]
        Km_Na = self._parameters["Km_Na"]
        PNaK = self._parameters["PNaK"]
        KmCa = self._parameters["KmCa"]
        KmNai = self._parameters["KmNai"]
        Ksat = self._parameters["Ksat"]
        alpha = self._parameters["alpha"]
        gamma = self._parameters["gamma"]
        kNaCa = self._parameters["kNaCa"]
        KPCa = self._parameters["KPCa"]
        g_PCa = self._parameters["g_PCa"]
        g_to = self._parameters["g_to"]
        Cao = self._parameters["Cao"]
        F = self._parameters["F"]
        Ki = self._parameters["Ki"]
        Ko = self._parameters["Ko"]
        Nao = self._parameters["Nao"]
        R = self._parameters["R"]
        T = self._parameters["T"]

        # Init return args
        current = [ufl.zero()]*1

        # Expressions for the Electric potentials component
        E_Na = R*T*ufl.ln(Nao/Nai)/F
        E_K = R*T*ufl.ln(Ko/Ki)/F
        E_Ks = R*T*ufl.ln((Ko + Nao*PkNa)/(Ki + PkNa*Nai))/F
        E_Ca = 0.5*R*T*ufl.ln(Cao/Cai)/F

        # Expressions for the i_Na component
        i_Na = g_Na*(m*m*m)*(-E_Na + V)*h*j

        # Expressions for the i_CaL component
        i_CaL = 4*g_CaL*ufl.elem_pow(F, 2.0)*(-0.341*Cao +\
            Cai*ufl.exp(2*F*V/(R*T)))*V*d*f1*f2*fCa/(R*T*(-1 +\
            ufl.exp(2*F*V/(R*T))))

        # Expressions for the i_Kr component
        i_Kr = 0.430331482912*g_Kr*ufl.sqrt(Ko)*(-E_K + V)*Xr1*Xr2

        # Expressions for the i_Ks component
        i_Ks = g_Ks*(Xs*Xs)*(1 + 0.6/(1 +\
            6.48182102606e-07*ufl.elem_pow(1.0/Cai, 1.4)))*(-E_Ks + V)

        # Expressions for the i_K1 component
        alpha_K1 = 3.91/(1 + 2.44592857399e-52*ufl.exp(594.2*V - 594.2*E_K))
        beta_K1 = (0.00277806676906*ufl.exp(588.6*V - 588.6*E_K) -\
            1.5394838221*ufl.exp(0.2*V - 0.2*E_K))/(1 + ufl.exp(454.7*V -\
            454.7*E_K))
        XK1_inf = alpha_K1/(alpha_K1 + beta_K1)
        i_K1 = 0.430331482912*g_K1*ufl.sqrt(Ko)*(-E_K + V)*XK1_inf

        # Expressions for the i_f component
        i_f = g_f*(-E_f + V)*Xf

        # Expressions for the i_b Na component
        i_b_Na = g_b_Na*(-E_Na + V)

        # Expressions for the i_b Ca component
        i_b_Ca = g_b_Ca*(-E_Ca + V)

        # Expressions for the i_NaK component
        i_NaK = Ko*PNaK*Nai/((Km_K + Ko)*(Km_Na + Nai)*(1 +\
            0.0353*ufl.exp(-F*V/(R*T)) + 0.1245*ufl.exp(-0.1*F*V/(R*T))))

        # Expressions for the i_NaCa component
        i_NaCa = kNaCa*(Cao*(Nai*Nai*Nai)*ufl.exp(F*gamma*V/(R*T)) -\
            alpha*(Nao*Nao*Nao)*Cai*ufl.exp(F*(-1 + gamma)*V/(R*T)))/((1 +\
            Ksat*ufl.exp(F*(-1 + gamma)*V/(R*T)))*(Cao +\
            KmCa)*((KmNai*KmNai*KmNai) + (Nao*Nao*Nao)))

        # Expressions for the i_PCa component
        i_PCa = g_PCa*Cai/(KPCa + Cai)

        # Expressions for the i_to component
        i_to = g_to*(-E_K + V)*q*r

        # Expressions for the Membrane component
        current[0] = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK -\
            i_PCa - i_b_Ca - i_b_Na - i_f - i_to

        # Return results
        return current[0]

    def I(self, v, s, time=None):
        """
        Transmembrane current

           I = -dV/dt

        """
        return -self._I(v, s, time)

    def F(self, v, s, time=None):
        """
        Right hand side for ODE system
        """
        time = time if time else Constant(0.0)

        # Assign states
        V = v
        assert(len(s) == 17)
        m, h, j, d, f1, f2, fCa, Xr1, Xr2, Xs, Xf, q, r, Nai, g, Cai, Ca_SR = s

        # Assign parameters
        g_Na = self._parameters["g_Na"]
        g_CaL = self._parameters["g_CaL"]
        tau_fCa = self._parameters["tau_fCa"]
        L0 = self._parameters["L0"]
        Q = self._parameters["Q"]
        g_b_Na = self._parameters["g_b_Na"]
        g_b_Ca = self._parameters["g_b_Ca"]
        Km_K = self._parameters["Km_K"]
        Km_Na = self._parameters["Km_Na"]
        PNaK = self._parameters["PNaK"]
        KmCa = self._parameters["KmCa"]
        KmNai = self._parameters["KmNai"]
        Ksat = self._parameters["Ksat"]
        alpha = self._parameters["alpha"]
        gamma = self._parameters["gamma"]
        kNaCa = self._parameters["kNaCa"]
        KPCa = self._parameters["KPCa"]
        g_PCa = self._parameters["g_PCa"]
        Cao = self._parameters["Cao"]
        Cm = self._parameters["Cm"]
        F = self._parameters["F"]
        Ko = self._parameters["Ko"]
        Nao = self._parameters["Nao"]
        R = self._parameters["R"]
        T = self._parameters["T"]
        V_SR = self._parameters["V_SR"]
        Vc = self._parameters["Vc"]
        Buf_C = self._parameters["Buf_C"]
        Buf_SR = self._parameters["Buf_SR"]
        Kbuf_C = self._parameters["Kbuf_C"]
        Kbuf_SR = self._parameters["Kbuf_SR"]
        Kup = self._parameters["Kup"]
        V_leak = self._parameters["V_leak"]
        VmaxUp = self._parameters["VmaxUp"]
        a_rel = self._parameters["a_rel"]
        b_rel = self._parameters["b_rel"]
        c_rel = self._parameters["c_rel"]
        tau_g = self._parameters["tau_g"]

        # Init return args
        F_expressions = [ufl.zero()]*17

        # Expressions for the Electric potentials component
        E_Na = R*T*ufl.ln(Nao/Nai)/F
        E_Ca = 0.5*R*T*ufl.ln(Cao/Cai)/F

        # Expressions for the i_Na component
        i_Na = g_Na*(m*m*m)*(-E_Na + V)*h*j

        # Expressions for the m gate component
        m_inf = ufl.elem_pow(1 + 0.00308976260789*ufl.exp(-169.491525424*V),\
            -0.333333333333)
        alpha_m = 1.0/(1 + 6.14421235333e-06*ufl.exp(-200.0*V))
        beta_m = 0.1/(1 + 1096.63315843*ufl.exp(200.0*V)) + 0.1/(1 +\
            ufl.exp(-1/4 + 5*V))
        tau_m = 0.001*alpha_m*beta_m
        F_expressions[0] = (-m + m_inf)/tau_m

        # Expressions for the h gate component
        h_inf = 1.0/ufl.sqrt(1 + 311490.091283*ufl.exp(175.438596491*V))
        alpha_h = ufl.conditional(ufl.lt(V, -0.04),\
            4.43126792958e-07*ufl.exp(-147.058823529*V), 0)
        beta_h = ufl.conditional(ufl.lt(V, -0.04), 310000.0*ufl.exp(348.5*V)\
            + 2.7*ufl.exp(79.0*V), 0.77/(0.13 +\
            0.0497581410839*ufl.exp(-90.0900900901*V)))
        tau_h = ufl.conditional(ufl.lt(V, -0.04), 1.5/(1000*alpha_h +\
            1000*beta_h), 0.002542)
        F_expressions[1] = (-h + h_inf)/tau_h

        # Expressions for the j gate component
        j_inf = 1/ufl.sqrt(1 + 311490.091283*ufl.exp(175.438596491*V))
        alpha_j = ufl.conditional(ufl.lt(V, -0.04), (37.78 +\
            1000*V)*(-25428.0*ufl.exp(244.4*V) -\
            6.948e-06*ufl.exp(-43.91*V))/(1 +\
            50262745826.0*ufl.exp(311.0*V)), 0)
        beta_j = ufl.conditional(ufl.lt(V, -0.04),\
            0.02424*ufl.exp(-10.52*V)/(1 +\
            0.0039608683399*ufl.exp(-137.8*V)), 0.6*ufl.exp(57.0*V)/(1 +\
            0.0407622039784*ufl.exp(-100.0*V)))
        tau_j = 7.0/(1000*alpha_j + 1000*beta_j)
        F_expressions[2] = (-j + j_inf)/tau_j

        # Expressions for the i_CaL component
        i_CaL = 4*g_CaL*ufl.elem_pow(F, 2.0)*(-0.341*Cao +\
            Cai*ufl.exp(2*F*V/(R*T)))*V*d*f1*f2*fCa/(R*T*(-1 +\
            ufl.exp(2*F*V/(R*T))))

        # Expressions for the d gate component
        d_infinity = 1.0/(1 + 0.42522244068*ufl.exp(-1000*V/7))
        alpha_d = 0.25 + 1.4/(1 + ufl.exp(-35/13 - 1000*V/13))
        beta_d = 1.4/(1 + ufl.exp(1 + 200*V))
        gamma_d = 1.0/(1 + ufl.exp(5/2 - 50*V))
        tau_d = 0.001*gamma_d + 0.001*alpha_d*beta_d
        F_expressions[3] = (-d + d_infinity)/tau_d

        # Expressions for the F1 gate component
        f1_inf = 1.0/(1 + 4485.77548844*ufl.exp(1000*V/3))
        constf1 = ufl.conditional(ufl.gt(-f1 + f1_inf, 0), 0.92835 +\
            1433.0*Cai, 1)
        tau_f1 = (20 + 180/(1 + 20.0855369232*ufl.exp(100.0*V)) + 200.0/(1 +\
            3.66929666762*ufl.exp(-100.0*V)) +\
            1102.5*ufl.exp(-0.00444444444444*ufl.elem_pow(27 + 1000*V,\
            4)))*constf1/1000
        F_expressions[4] = (-f1 + f1_inf)/tau_f1

        # Expressions for the F2 gate component
        f2_inf = 0.33 + 0.67/(1 + 2456.51756009*ufl.exp(250*V))
        constf2 = 2.0
        tau_f2 = constf2*(16/(1 + 20.0855369232*ufl.exp(100.0*V)) + 31/(1 +\
            12.1824939607*ufl.exp(-100.0*V)) +\
            600*ufl.exp(-0.00588235294118*((25 + 1000*V)*(25 +\
            1000*V))))/1000
        F_expressions[5] = (-f2 + f2_inf)/tau_f2

        # Expressions for the FCa gate component
        alpha_fCa = 1.0/(1 + 5.95374180765e+25*ufl.elem_pow(Cai, 8))
        beta_fCa = 0.1/(1 + 0.000123409804087*ufl.exp(10000.0*Cai))
        gamma_fCa = 0.3/(1 + 0.391605626677*ufl.exp(1250.0*Cai))
        fCa_inf = 0.760109455762*alpha_fCa + 0.760109455762*beta_fCa +\
            0.760109455762*gamma_fCa
        constfCa = ufl.conditional(ufl.And(ufl.gt(V, -0.06), ufl.gt(fCa_inf,\
            fCa)), 0, 1)
        F_expressions[6] = (-fCa + fCa_inf)*constfCa/tau_fCa

        # Expressions for the Xr1 gate component
        V_half = -19.0 - 1000*R*T*ufl.ln(ufl.elem_pow(1 + 0.384615384615*Cao,\
            4)/(L0*ufl.elem_pow(1 + 1.72413793103*Cao, 4)))/(F*Q)
        Xr1_inf = 1.0/(1 + ufl.exp(0.204081632653*V_half - 204.081632653*V))
        alpha_Xr1 = 450.0/(1 + ufl.exp(-9/2 - 100*V))
        beta_Xr1 = 6.0/(1 + 13.5813245226*ufl.exp(86.9565217391*V))
        tau_Xr1 = 0.001*alpha_Xr1*beta_Xr1
        F_expressions[7] = (-Xr1 + Xr1_inf)/tau_Xr1

        # Expressions for the Xr2 gate component
        Xr2_infinity = 1.0/(1 + ufl.exp(44/25 + 20*V))
        alpha_Xr2 = 3.0/(1 + ufl.exp(-3 - 50*V))
        beta_Xr2 = 1.12/(1 + ufl.exp(-3 + 50*V))
        tau_Xr2 = 0.001*alpha_Xr2*beta_Xr2
        F_expressions[8] = (-Xr2 + Xr2_infinity)/tau_Xr2

        # Expressions for the Xs gate component
        Xs_infinity = 1.0/(1 + ufl.exp(-5/4 - 125*V/2))
        alpha_Xs = 1100.0/ufl.sqrt(1 + ufl.exp(-5/3 - 500*V/3))
        beta_Xs = 1.0/(1 + ufl.exp(-3 + 50*V))
        tau_Xs = 0.001*alpha_Xs*beta_Xs
        F_expressions[9] = (-Xs + Xs_infinity)/tau_Xs

        # Expressions for the Xf gate component
        Xf_infinity = 1.0/(1 + 5780495.71031*ufl.exp(200*V))
        tau_Xf = 1.9/(1 + ufl.exp(3/2 + 100*V))
        F_expressions[10] = (-Xf + Xf_infinity)/tau_Xf

        # Expressions for the i_b Na component
        i_b_Na = g_b_Na*(-E_Na + V)

        # Expressions for the i_b Ca component
        i_b_Ca = g_b_Ca*(-E_Ca + V)

        # Expressions for the i_NaK component
        i_NaK = Ko*PNaK*Nai/((Km_K + Ko)*(Km_Na + Nai)*(1 +\
            0.0353*ufl.exp(-F*V/(R*T)) + 0.1245*ufl.exp(-0.1*F*V/(R*T))))

        # Expressions for the i_NaCa component
        i_NaCa = kNaCa*(Cao*(Nai*Nai*Nai)*ufl.exp(F*gamma*V/(R*T)) -\
            alpha*(Nao*Nao*Nao)*Cai*ufl.exp(F*(-1 + gamma)*V/(R*T)))/((1 +\
            Ksat*ufl.exp(F*(-1 + gamma)*V/(R*T)))*(Cao +\
            KmCa)*((KmNai*KmNai*KmNai) + (Nao*Nao*Nao)))

        # Expressions for the i_PCa component
        i_PCa = g_PCa*Cai/(KPCa + Cai)

        # Expressions for the q gate component
        q_inf = 1.0/(1 + ufl.exp(53/13 + 1000*V/13))
        tau_q = 0.00606 + 0.039102/(0.0168716780457*ufl.exp(-80.0*V) +\
            6.42137321286*ufl.exp(100.0*V))
        F_expressions[11] = (-q + q_inf)/tau_q

        # Expressions for the r gate component
        r_inf = 1.0/(1 + 3.28489055021*ufl.exp(-53.3333333333*V))
        tau_r = 0.00275352 + 0.01440516/(16.3010892258*ufl.exp(90.0*V) +\
            0.0211152735604*ufl.exp(-120.0*V))
        F_expressions[12] = (-r + r_inf)/tau_r

        # Expressions for the Sodium dynamics component
        F_expressions[13] = -1e+18*Cm*(3*i_NaCa + 3*i_NaK + i_Na +\
            i_b_Na)/(F*Vc)

        # Expressions for the Calcium dynamics component
        i_rel = 0.0556*(c_rel + a_rel*(Ca_SR*Ca_SR)/((b_rel*b_rel) +\
            (Ca_SR*Ca_SR)))*d*g
        i_up = VmaxUp/(1 + (Kup*Kup)/(Cai*Cai))
        i_leak = V_leak*(-Cai + Ca_SR)
        g_inf = ufl.conditional(ufl.le(Cai, 0.00035), 1.0/(1 +\
            5.43991024148e+20*ufl.elem_pow(Cai, 6)), 1.0/(1 +\
            1.9720198874e+55*ufl.elem_pow(Cai, 16)))
        const2 = ufl.conditional(ufl.And(ufl.gt(V, -0.06), ufl.gt(g_inf, g)),\
            0, 1)
        F_expressions[14] = (-g + g_inf)*const2/tau_g
        Cai_bufc = 1.0/(1 + Buf_C*Kbuf_C/((Kbuf_C + Cai)*(Kbuf_C + Cai)))
        Ca_SR_bufSR = 1.0/(1 + Buf_SR*Kbuf_SR/((Kbuf_SR + Ca_SR)*(Kbuf_SR +\
            Ca_SR)))
        F_expressions[15] = (-i_up - 5e+17*Cm*(-2*i_NaCa + i_CaL + i_PCa +\
            i_b_Ca)/(F*Vc) + i_leak + i_rel)*Cai_bufc
        F_expressions[16] = Vc*(-i_leak - i_rel + i_up)*Ca_SR_bufSR/V_SR

        # Return results
        return dolfin.as_vector(F_expressions)

    def num_states(self):
        return 17

    def __str__(self):
        return 'Paci2013_atrial cardiac cell model'
