from math import *
import matplotlib.pyplot as plt


g = 9.81+1.2
SafetyFact = 2.1 ##for average materials operated in ordinary environments and sub-jected to loads and stresses that can be determined.
Sy = 450000000
Su = 600000000
Ppa = 3700

Ma = 0
Mb = 0
Mc = 0
Md = 0
###################################################################################### chemin clutch du bas fermé
def engren1a():
    Rp = 0.05
    phi = 0.349 ##alpha engrenage normal,20°
    gamma = 0.78 ##45°
    Re = 0.15
    Lab = 0.3
    Lbc = 0.3
    Lcd = 0.3
    Ma = 0
    Mc = 0
    poidpassager = Ppa*1 #kg
    Pasc = poidpassager * g
    Pcont = poidpassager/2*g
    Fyb = -(Pasc + Pcont)
    Mb = -Rp*g*poidpassager/2
    print(Mb)
    Md = -Mb
    Fxa = 0
    Fxd = 0
    Fxb = 0 ##Hypothèse
    Fxc = 0 ##Hypothèse
    Ft = 0
    Fyd = 0
    Fyc = (-Fyb*Lab - Fyd*(Lab+Lbc+Lcd))/(Lbc+Lab)
    Fya = -Fyb-Fyc-Fyd
    Fzd = 0
    Fzb = 0 ##hypothèse
    Fzc = (-Fzb*Lab - Fzd*(Lab+Lbc+Lcd))/(Lbc+Lab)
    Fza = -Fzb-Fzc-Fzd
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Ma,Mc,Md,Mb,Rp,phi,gamma,Re,Lab,Lbc,Lcd

#Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Ma,Mc,Md,Mb,Rp,phi,gamma,Re,Lab,Lbc,Lcd = engren1()




def engren2a(M):
    Rp = 0.15
    phi = 0.349 ##alpha engrenage normal,20°
    gamma = 0.78 ##45°
    Re = 0.1
    Lab = 0.3
    Lbc = 0.3
    Lcd = 0.3
    Ma = -M
    Mc = 0
    Mb = 0
    Md = -Ma
    Fxb = Md/Re*tan(phi)*sin(gamma)
    Fxd = -Fxb
    Fxa = 0 ##Hypothèse
    Fxc = 0 ##Hypothèse
    Ft = Md/Re
    Fya = 0
    Fyd = abs(Ft)*tan(phi)*cos(gamma)
    Fyc = (Fya*Lab - Fyd*(Lbc+Lcd) - Fxd*Re)/Lbc
    Fyb = -Fya-Fyc-Fyd
    Fzd = Ft
    Fza = 0
    Fzb = 0 ##hypothèse
    Fzc = (Fza*Lab-Fzd*(Lbc+Lcd))/Lbc
    Fzb = -Fza-Fzc-Fzd
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Ma,Mc,Md,Mb,Rp,phi,gamma,Re,Lab,Lbc,Lcd

#Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Ma,Mc,Md,Mb,Rp,phi,gamma,Re,Lab,Lbc,Lcd = engren2a(Md)


###########################################################chemin clutch du haut fermé

def engren1c():
    alpha = 0.78 ##45°
    Rp = 0.05
    Re = 0.15
    Lab = 0.3
    Lbc = 0.3
    Lcd = 0.3
    Ma = 0
    Md = 0
    poidpassager = Ppa #kg
    Pasc = poidpassager * g
    Pcont = poidpassager/2*g
    Fyb = -(Pasc + Pcont)
    Mb = -Rp*g*poidpassager/2
    Mc = -Mb
    Fxa = 0 ##Hypothèse
    Fxd = 0 ##Hypothèse
    Fxb = 0 ##Hypothèse
    Fxc = 0 ##Hypothèse
    Ft = Mc/Re
    Fyc = -abs(Ft)*tan(alpha)
    Fya = (-Fyb*(Lbc+Lcd) - Fyc*Lcd)/(Lab+Lbc+Lcd)
    Fyd = -Fya-Fyc-Fyb
    Fzc = -Ft
    print(Mc)
    Fzb = 0 ##hypothèse
    Fza = (-Fzb*(Lbc+Lcd) - Fzc*Lcd)/(Lab+Lbc+Lcd)
    Fzd = -Fza-Fzb-Fzc
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Ma,Mc,Md,Mb,Rp,Re,Lab,Lbc,Lcd

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Ma,Mc,Md,Mb,Rp,Re,Lab,Lbc,Lcd = engren1c()

##mettre Mc dedans
def engren3c(M):
    alpha = 0.78 ##45°
    Re = 0.15
    Lab = 0.3
    Lbc = 0.3
    Mb2 = M
    Fyb = 0
    Mb = 0
    Mc = 0
    Ma = 0
    Fxa = 0 ##Hypothèse
    Fxb = 0 ##Hypothèse
    Fxc = 0 ##Hypothèse
    Fxd = 0
    Fyc = 0
    Fya = 0
    Fyb = 0
    Fyd = 0
    Md = 0
    print(M)
    Fzb = 2*M/Re
    Fzc = -Fzb*Lab/(Lab+Lbc)
    Fza = -Fzb-Fzc
    Fzd = 0
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ma,Mc,Mb,Md,Rp,alpha,Re,Lab,Lbc,Lcd

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ma,Mc,Mb,Md,Rp,alpha,Re,Lab,Lbc,Lcd = engren3c(Mc)

##mettre Fzb/2
def engren4c(F):
    Re2 = 0.15
    alpha1 = 0.349 ##alpha engrenage normal,20°
    Lab = 0.3
    Lbc = 0.3
    Lcd = 0.3
    Mc = 0
    Ma = 0
    Mb = -F*Re2
    Md = -Mb
    print(Md)
    Fzb = -F
    Fyb = tan(alpha1)*abs(F)
    Fzd = 0
    Fxa = 0 ##hypothèse
    Fyd = 0
    Fxb = 0 ##hypothèse
    Fxc = 0 ##hypothèse
    Fxd = 0 ##hypothèse
    Fyc = (-Fyb*Lab - Fyd*(Lab+Lbc+Lcd))/(Lbc+Lab)
    Fya = -Fyb-Fyc-Fyd
    Fzc = (-Fzb*Lab - Fzd*(Lab+Lbc+Lcd))/(Lbc+Lab)
    Fza = -Fzb-Fzc-Fzd
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Md,Mb,Mc,Ma,Re2,alpha1,Lab,Lbc,Lcd

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Md,Mb,Mc,Ma,Re2,alpha1,Lab,Lbc,Lcd = engren4c(Fzb/2)

##mettre Md
def engren5c(M):
    Re2 = 0.45/2
    alpha1 = 0.349 ##alpha engrenage normal,20°
    Lab = 0.3
    Lbc = 0.3
    Lcd = 0.3
    Ma = -M
    Mb = 0
    Md = 0
    Mc = -Ma
    Fza = 0
    Fxa = 0 ##hypothèse
    Fya = 0
    Fxb = 0 ##hypothèse
    Fxc = 0 ##hypothèse
    Fxd = 0 ##hypothèse
    Fzc = Mc/Re
    Fyc = tan(alpha1)*abs(Fzc)
    Fyb = (-Fya*(Lab+Lbc+Lcd) - Fyc*Lcd)/(Lbc+Lcd)
    Fyd = -Fya-Fyc-Fyb
    Fzd = (Fza*Lab - Fzc*Lbc)/(Lbc+Lcd)
    Fzb = -Fza-Fzc-Fzd
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Md,Mc,Mb,Ma,Re2,alpha1,Lab,Lbc,Lcd

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Md,Mc,Mb,Ma,Re2,alpha1,Lab,Lbc,Lcd = engren5c(Md)

def engren2c(M):
    print(M)
    Re2 = 0.45/2
    Re3 = 0.1
    phi = 0.349 ##alpha engrenage normal,20°
    alpa = 0.349 ##alpha engrenage normal,20°
    gamma = 0.78 ##45°
    Lab = 0.3
    Lbc = 0.3
    Lcd = 0.3
    Ma = 0
    Mb = M
    Md = -Mb
    Mc = 0
    Fzb = -Mb/Re2
    Fzd = Md/Re3
    print(Fzd)
    Fzc = (-Fzb*Lab - Fzd*(Lab+Lbc+Lcd))/(Lbc+Lab)
    Fza = -Fzb-Fzc-Fzd
    Fxa = abs(Fzd)*tan(phi)*sin(gamma)
    Fxb = 0 ##hypothèse
    Fxc = 0 ##hypothèse
    Fxd = -Fxa ##le support de force axial est en a
    Fyd = tan(phi)*cos(gamma)*abs(Fzd)
    Fyb = tan(alpha)*Fzb
    Fyc = (-Fyb*Lab - Fyd*(Lab+Lbc+Lcd) - Fxd*Re3)/(Lbc+Lab)
    Fya = -Fyd-Fyc-Fyb
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Md,Mc,Mb,Ma,Re2,alpha,Lab,Lbc,Lcd

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Md,Mc,Mb,Ma,Re2,alpha,Lab,Lbc,Lcd = engren2c(Mc)




#############################################################################

##Ft et Ft2, force tangente venant des deux engrenages, peuvent avoir des signes différents en fonction de la configuration
##positif si force tangente vers le fond
def engren6(Ft,Ft2=0):
    Re2 = 0.1
    Re3 = 0.2
    phi = 0.349 ##alpha engrenage normal,20°
    gamma = 0.78 ##45°
    alpha = 0.349 ##alpha engrenage normal,20°
    Lab = 0.2
    Lbc = 0.2
    Lcd = 0.2
    Ma = (-Ft+Ft2)*Re2
    Md = -Ma
    Mc = 0
    Mb = 0
    Fya = (abs(Ft)-abs(Ft2))*tan(phi)*cos(gamma)
    Fzd = -Md/Re3
    Fza = -Ft-Ft2
    Fxa1 = abs(Ft)*tan(phi)*sin(gamma)
    Fxa2 = abs(Ft2)*tan(phi)*sin(gamma)
    Fxa = Fxa1 + Fxa2
    Fyd = -tan(alpha)*abs(Fzd)
    Fxb = 0 ##hypothèse
    Fxc = -Fxa
    Fxd = 0 ##hypothèse
    Fyc = (Fya*Lab - Fyd*(Lbc+Lcd) -(Fxa1-Fxa2)*Re2)/Lbc
    Fyb = -Fya-Fyc-Fyd
    Fzc = (Fza*Lab-Fzd*(Lbc+Lcd))/Lbc
    Fzb = -Fza-Fzc-Fzd
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Md,Ma,Mc,Mb,Re2,phi,gamma,Re3,Lab,Lbc,Lcd,alpha

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Md,Ma,Mc,Mb,Re2,phi,gamma,Re3,Lab,Lbc,Lcd,alpha = engren6(Fzd)
##mettre Fzd en entrée

##partie moteur
def engren7(F):
    Re2 = 0.04
    alpha1 = 0.349 ##alpha engrenage normal,20°
    Lab = 0.2
    Lbc = 0.2
    Lcd = 0.2
    Ma = -F*Re2
    Md = -Ma
    Fza = -F
    Fya = tan(alpha1)*abs(F)
    Fzd = 0 ##hypothèse
    Fxa = 0 ##hypothèse
    Fyd = 0 ##hypothèse
    Fxb = 0 ##hypothèse
    Fxc = 0 ##hypothèse
    Fxd = 0 ##hypothèse
    Fyc = (Fya*Lab - Fyd*(Lbc+Lcd))/Lbc
    Fyb = -Fya-Fyc-Fyd
    Fzc = (Fza*Lab-Fzd*(Lbc+Lcd))/Lbc
    Fzb = -Fza-Fzc-Fzd
    return Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Md,Ma,Re2,alpha1,Lab,Lbc,Lcd

Fya,Fyb,Fyc,Fyd,Fxa,Fxb,Fxc,Fxd,Fza,Fzb,Fzc,Fzd,Ft,Md,Ma,Re2,alpha1,Lab,Lbc,Lcd = engren7(Fzd)




#######graphe et calcul de diamètres

Xstart = 0
Xend = Lab + Lbc + Lcd
nbrstep = 500
dX = (Xend - Xstart)/nbrstep
x = 0
X = []
Fx = []
Fy = []
Fz = []
T = [] ## la torsion est positive si elle tourne le shaft vers le contre poids
My = [] ##le couple positif est celui qui étire les fibres du coté inférieur
Mz = [] ##le couple positif est celui qui étire les fibres du coté du contrepoid
Ftot = []
Mtot = []
Tmax = -10
Ftotmax = -10
Mtotmax = -10
Fxmax = -10


for i in range(0,nbrstep):
    X.append(x)
    fx = 0
    fy = 0
    fz = 0
    t = 0
    ftot = 0
    mtot = 0
    my = 0
    mz = 0
    if x < Lab:
        fy = -Fya
        fx = -Fxa
        fz = -Fza
        my = x*Fya
        mz = x*Fza
        t = -Ma
        
    elif Lab <= x and x < Lab + Lbc:
        fy = -Fya -Fyb
        fx = -Fxa -Fxb
        fz = -Fza -Fzb
        my = x*Fya + (x-Lab)*Fyb
        mz = x*Fza + (x-Lab)*Fzb
        t = -Ma-Mb
        
    elif Lab + Lbc <= x and x <= Lab + Lbc + Lcd:
        fy = -Fya -Fyb - Fyc
        fx = -Fxa -Fxb - Fxc
        fz = -Fza -Fzb - Fzc
        my = x*Fya + (x-Lab)*Fyb + (x-(Lab+Lbc))*Fyc
        mz = x*Fza + (x-Lab)*Fzb + (x-(Lab+Lbc))*Fzc
        t = -Ma-Mb-Mc
    
        
    x += dX
    ftot = sqrt((fx**2)+(fy**2)+(fz**2))
    mtot = sqrt((my**2)+(mz**2))
    if Ftotmax < ftot :
        Ftotmax = ftot
    if Mtotmax < mtot :
        Mtotmax = mtot
    if Tmax < abs(t) :
        Tmax = t
    if Fxmax < abs(fx) :
        Fxmax = abs(fx)
    
    Fx.append(fx)
    Fy.append(fy)
    Fz.append(fz)
    T.append(t)
    My.append(my)
    Mz.append(mz)
    Ftot.append(ftot)
    Mtot.append(mtot)
    
print("debut")
print(Fxa)
print(Fxb)
print(Fxc)
print(Fxd)
print(sqrt(Fza**2 + Fya**2))
print(sqrt(Fzb**2 + Fyb**2))
print(sqrt(Fzc**2 + Fyc**2))
print(sqrt(Fzd**2 + Fyd**2))
plt.title("forces en fonction de la distance", fontsize=8)
plt.xlabel("distance [m]", fontsize=8)
plt.ylabel("force [N]", fontsize=8)
plt.plot(X,Fx,label = "Fx")
plt.plot(X,Fy,label = "Fy")
plt.plot(X,Fz,label = "Fz")
plt.plot(X,My,label = "My")
plt.plot(X,Mz,label = "Mz")
plt.plot(X,T,label = "T")
plt.plot(X,Ftot,label = "Ftot")
#plt.plot(X,Mtot,label = "Mtot")
plt.legend(loc="upper left", fontsize=8)
plt.grid()
print("valeur de force max :"+ str(Ftotmax) +"N")
print("valeur de couple max :"+ str(Mtotmax) +"N/m")
print("valeur de torsion max :"+ str(Tmax) +"N/m")
#plt.savefig("exo5-4.png")
plt.show()
plt.clf()

dmin = ((((32*Mtotmax)/pi)**2 + 3*((16*Tmax/pi)**2))*((SafetyFact/Sy)**2))**(1/6)
print("valeur de dmin :"+ str(dmin) +"m")


###fatigue
guessd = 0.001
deltad = 0.0001
CL = 1
CG = 0.9
CS = 0.77
CT = 1
CR = 0.753

Su = 600000000 #MPa
Sn= 0.5*Su*CL*CG*CS*CT*CR
Kf_sigma = 1 + 0.725*(1.3-1)
Kf_tau = 1 + 0.76*(1.12-1)

Fa = Fxmax
d = guessd
print(Mtotmax)
print(Tmax)
print(Fa)
print(Sn)
sigmaEA = 32*Mtotmax/(pi*d**3)*Kf_sigma
sigmaEM = 8*Fa/((pi**2) * d**4)*Kf_sigma + sqrt((16*Tmax/(pi*(d**3))*Kf_tau)**2 + (8*Fa/(pi**2 * d**4)*Kf_sigma)**2 )
while sigmaEA > -(Sn/Su)*sigmaEM + Sn or sigmaEA > - sigmaEM + Sy:
    d += deltad
    sigmaEA = 32*Mtotmax/(pi*d**3)*Kf_sigma
    sigmaEM = 8*Fa/(pi**2 * d**4)*Kf_sigma + sqrt((16*Tmax/(pi*d**3)*Kf_tau)**2 + (8*Fa/(pi**2 * d**4)*Kf_sigma)**2 )
    
#rapport = 2*Mtotmax*Kf_sigma/(Tmax*Kf_tau)
#print(rapport)
#sigma_m_max = Sn/(Sn/Su + rapport)

#d = (16*Tmax*Kf_tau/(pi*sigma_m_max))**(1/3)
print("valeur de d avec fatigue:"+ str(d) +"m")
