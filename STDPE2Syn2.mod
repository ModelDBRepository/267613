: STDP by Hines, changed to dual exponential (BPG 6-1-09)
: Modified by BPG 13-12-08
: Limited weights: max weight is wmax and min weight is wmin
: (initial weight is specified by netconn - usually set to wmin)
: Rhythmic GABAB suppresses conductance and promotes plasticity.
: When GABAB is low, conductance is high and plasticity is off.

NEURON {
	POINT_PROCESS STDPE2bis
	RANGE tau1, tau2, e, i, d, p, dtau, ptau, thresh, wmax, wmin, srcnt1, srcnt2, srcnt3, srcnt4
	RANGE g, gbdel, gblen, gbint, gscale, factor,dshift,dM,dV,B,C
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	srcnt1=0
	srcnt2=0
	srcnt3=0
	srcnt4=0
	tau1=.5 (ms) <1e-9,1e9>
	tau2 = 3 (ms) <1e-9,1e9>
	e = 0	(mV)
     pi=3.14159
	wmax = 0.008(uS)
      :wmax = 0.004 (uS)
	wmin = 0(uS)	: not used - use netconn weight instead (BPG)
	d = 0.3 :depression factor (multiplicative to prevent < 0) 0.8
	p = 1 : potentiation factor (additive, non-saturating) 1
	
       dM = -24   (ms) :-22
       dV= 6.32 :5    (ms)
	 ptau = 2 (ms) : Nishiyama2000 10

	thresh = -20 (mV)	: postsynaptic voltage threshold
     : thresh = -55 (mV)	: postsynaptic voltage threshold
	gbdel = 50 (ms) <1e-9,1e9> : initial GABAB off interval (ms)
	gbint = 125 (ms) <1e-9,1e9> : GABAB off interval (ms)
	gblen = 125 (ms) <1e-9,1e9> : GABAB on length (ms)
	gscale = 1	: relative suppression by GABAB
      dshift=0 (ms)
}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
	on
	g (uS)
	gs
	factor
}

STATE {
	C (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	C = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	gs=1
	on=1	: initially not plastic
	tpost =-1e9
	net_send(0, 1)
	:net_send(gbdel, 3)	: initial GABAB off period
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - C
     	i = g*gs*(v - e)
    
}

DERIVATIVE state {
	C' = -C/tau1
	B' = -B/tau2
}

NET_RECEIVE(w (uS), A, tpre (ms) ) {
	INITIAL { A = 0.0  tpre = -1e9 }
	if (flag == 0) { : presynaptic spike  (after last post so depress)
		:printf("presynaptic spike  t=%g w=%g A=%g tpre=%g tpost=%g\n", t, w, A, tpre, tpost)
:		g = g + w + A	: only for single exp (BPG)
		C = C + (w + A)*factor
		B = B + (w + A)*factor
          	:printf(" B %f\t C %f\t w %f\n",B, C, w)
		tpre = t
		:if (srcnt1 !=0) {printf("presynaptic_spike apical_%g(%g) syn(%g)(%g) tpost= %g tpre= %g\n", srcnt1, srcnt2, srcnt4, srcnt3, tpost, t)}
		if (on == 1) {
			A = A * (1-(d*exp(-((tpost-t)-dM)^2/(2*dV*dV))))
			:printf("%g %g\n", t, (1-(d*exp(-((tpost-t)-dM)^2/(2*dV*dV)))))
		}
	}else if (flag == 2 && on == 1) { : postsynaptic spike
		:if (srcnt1 !=0) {printf("postsynaptic spike apical_%g(%g) t=%g tpost=%g\n", srcnt1, srcnt2, t, tpost)}
		tpost = t
            FOR_NETCONS(w1, A1, tp) { : also can hide NET_RECEIVE args
			:printf("entry FOR_NETCONS w1=%g A1=%g tp=%g t=%g\n", w1, A1, tp, t)
			:if (srcnt1 !=0) {printf("postsynaptic_spike apical_%g(%g) syn(%g)(%g) tpost= %g tpre= %g\n", srcnt1, srcnt2, srcnt4, srcnt3, t, tp)}
			A1 = A1 + (2*w1-w1-A1)*p*exp((tp - t)/ptau)
			:printf("%g %g\n", t, A1)
			}
	} else if (flag == 1) { : flag == 1 from INITIAL block
		:printf("entry flag=%g t=%g\n", flag, t)
		WATCH (v > thresh) 2
	}
	else if (flag == 3) { : plasticity control
		if (on == 0) { : start plasticity
			on = 1
			gs = gscale
			net_send(gblen, 3)
		}
		else { : end burst
			on = 0
			gs = 1
			net_send(gbint, 3)
		}
	}
}

