import math
from numpy import log as ln
#####IEC60949-1988 
##Calculation of thermal permissible short-circuit current
## taking into account non-adiabatic heating effects############


####Short circuit current, and temperature
def I_short(F_e, I_AD):
	## F_e : factor to allow for heat loss in to the adjacent components
	## I_AD : short circuit current calculate at adiabatic basis
	## I_short : permissible short-circuit current
	I_short=F_e*I_AD
	return I_short

def I_AD(t, K, CSA, B, T_i, T_f):
	## t : duration
	## K : constant depending on the material of current carryin component (IEC Table I)
	## SCA : geometrical CAS of conductor
	## B : reciprocal of temperature coeffient of resistance of current carrying component at 0degree (table I)
	## T_i, T_f: initial and final temperature
	A=(K**2)*(CSA**2)*ln((T_f+B)/(T_i+B))
	I_AD=((A/t)**0.5)/1000
	return I_AD

def T_short(T_i,B,K,I_AD,CSA,t):
	A=((I_AD*1000)**2)*t
	C=(K**2)*(CSA**2)
	D=A/C
	T_short=(T_i+B)*math.exp(D)-B
	return T_short

### Factor to allow for heat loss in to the adjacent components

def F_e(t,CSA,F,SpH_C,SpH_insu,Rth_in):
	## X,Y : incorporating the thermal contact factor of 0.7 (1.0 for oil filled cable) See table III
	## SpH_C, SpH_insu, specifi heat of conductor and the adjacent insulator
	## Rth_: thermal resistivity of adjacent insulation
	A=(2464/SpH_C)*((SpH_insu/Rth_in)**0.5)
	B=(1.22//SpH_C)*(SpH_insu/Rth_in)
	D=t/CSA
	E=D**0.5
	C=1+F*A*E+(F**2)*B*D
	F_e=C**0.5
	return F_e

### Factor to allow for sheaths, screen and wires

def M(SpH_1, SpH_2,SpH_3,Rth_2,Rth_3,Thickness,F):
	## M : Sheath thermal contact factor
	## SpH_1 : Specific heat of sheath or screen
	## SpH_2, SpH_3: Specific heat of adjacent material
	## Rth_2,Rth_3: Thermal resistivity of on either side of sheath
	## Thickness: sheath thichness
	## F : thermal contact impact factor, 0.7 in general, but 0.9 when 
	     ## metalic component is completely bonded on one side to the adjacent medium
	A= (SpH_2/Rth_2)**0.5
	B= (SpH_3/Rth_3)**0.5
	C= 2*SpH_1*Thickness/1000
	M=(A+B)*F/C
	return M
def F_e_sheath(M, t):
	t_root=t**0.5
	A=0.61*M*t_root
	B=0.069*(A**2)
	C=0.0043*(A**3)
	F_e_sheath=1+A+B+C
	return F_e_sheath