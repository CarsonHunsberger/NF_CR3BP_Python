import numpy as np
from scipy.io import loadmat
from scipy.integrate import solve_ivp

N=11

datastr = 'NF_CR3BP/data/'

# mu = 0.012154
# mustr = '012154'

mu = 0.01215058560962404068939157752993196481839
mustr = '012150585609624041'

# mu = 0.0000030542
# mustr = '0541999999999999e-06'




Nstr = 'N11'
nftypes = ['Birkhoff','Resonant']
lpts = ['L1','L2','L3']
NFdata = {}
for n in range(len(nftypes)):
    typedata = {}
    for k in range(len(lpts)):
        tempstr = datastr+mustr+'/'+Nstr+'/'+nftypes[n]+'/'+lpts[k]+'/'
        anlqpNtoqp0 = loadmat(tempstr+'anlqpNtoqp0.mat')
        anlqpNtoqp0 = anlqpNtoqp0['anlqpNtoqp0']

        anlqp0toqpN = loadmat(tempstr+'anlqp0toqpN.mat')
        anlqp0toqpN = anlqp0toqpN['anlqp0toqpN']

        HAA = loadmat(tempstr+'HAA.mat')
        HAA = HAA['HAA']

        GenFuncEOMs = loadmat(tempstr+'GenFuncEOMs.mat')
        GenFuncEOMs = GenFuncEOMs['GenFuncEOMs']

        gamma = loadmat(tempstr+'gamma.mat')
        gamma = gamma['gamma']
        gamma = gamma[0,0]

        if k==0:
            aLpt = 1-gamma
            gammascale = 1
        if k==1:
            aLpt = 1+gamma
            gammascale = 1
        if k==2:
            aLpt = -gamma
            gammascale = -1

        
        C = loadmat(tempstr+'C.mat')       
        C = C['C']
        Cinv = np.linalg.inv(C)

        T1 = gammascale*gamma*np.eye(6,6)
        T1[2,2] = gamma
        T1[5,5] = gamma         
        T1inv = (1/gamma)*np.eye(6,6)

        temparr = {'anlqpNtoqp0':anlqpNtoqp0,'anlqp0toqpN':anlqp0toqpN,'GenFuncEOMs':GenFuncEOMs,
                   'gamma':gamma,'aLpt':aLpt,'gammascale':gammascale,'C':C,'Cinv':Cinv,'T1':T1,'T1inv':T1inv,'HAA':HAA}
        
        AApartialscell = loadmat(tempstr+'AApartialscell.mat')
        AApartialscell = AApartialscell['AApartialscell']
        temparr['AApartialscell'] = AApartialscell

        typedata[lpts[k]] = temparr
    NFdata[nftypes[n]] = typedata

optsnum = {'rtol': 1e-12,'atol': 1e-16}
opts = {'rtol': 1e-13,'atol': 1e-22}

temparr1 = np.linspace(2,N-1,N-2,dtype=int)
temparr = np.flip(temparr1)


V = np.eye(6,6)                 #Same for everything, can keep this here
V[4,0] = -1
V[3,1] = 1
Vinv = np.linalg.inv(V)



def celleval(x,cell):
   if len(cell[0]):
        return np.sum(cell[0].T * np.prod((x ** cell[1]),axis=1),axis=1)
   else:
      return np.zeros(6)

def AAtotilde(x,nftype):
    tilde = np.zeros(6)
    if nftype=='Resonant':
        x = x.copy()
        x[2] = x[2]-x[1]
        x[4] = x[4]+x[5]
    tilde[0] = np.sqrt(x[0].real)*np.exp(x[3].real)
    tilde[1] = np.sqrt(2*x[1].real)*np.cos(x[4].real)
    tilde[2] = np.sqrt(2*x[2].real)*np.cos(x[5].real)
    tilde[3] = np.sqrt(x[0].real)*np.exp(-x[3].real)
    tilde[4] = -np.sqrt(2*x[1].real)*np.sin(x[4].real)
    tilde[5] = -np.sqrt(2*x[2].real)*np.sin(x[5].real)
    # if (x[3].imag == 0):
    if (x[3].imag == 0.5*np.pi):
        tilde[0] = -tilde[0]
    elif (x[3].imag == np.pi):
        tilde[0] = -tilde[0]
        tilde[3] = -tilde[3]
    else:
        tilde[3] = -tilde[3]
    return tilde.real
   
def tildetoAA(x,nftype):
    AA = np.zeros(6,complex)
    AA[0] = np.abs(x[0]*x[3])
    AA[1] = 0.5*(np.power(x[1],2)+np.power(x[4],2))
    AA[2] = 0.5*(np.power(x[2],2)+np.power(x[5],2))
    if x[3] == 0:
        AA[3] = 0
    else:
        AA[3] = np.log(np.abs(x[0]/x[3]))

    if x[0] < 0:
        if x[3] < 0:
            AA[3] = AA[3]+(1j*np.pi)
        else:
            AA[3] = AA[3]+(0.5j*np.pi)
    else:
        if x[3] < 0:
            AA[3] = AA[3]+(1.5j*np.pi)
    AA[4] = -np.arctan2(x[4],x[1])
    AA[5] = -np.arctan2(x[5],x[2])
    if nftype=='Resonant':
        AA[2] = AA[1] + AA[2]
        AA[4] = AA[4] - AA[5]
    if AA[4].real < 0:
        AA[4] = 2*np.pi+AA[4].real
    if AA[5].real < 0:
        AA[5] = 2*np.pi+AA[5].real

    return AA

def sixdeval(x,cells):
    result = np.zeros(6,complex)
    for n in range(6):
        cell = cells[0,n][0]
        temp = celleval(x,cell)
        result[n] = temp[0]
    return result.real

def resAApartials(AA,cells):
    result = np.zeros(6)
    x = np.zeros(3)
    AA = AA.real
    x[0],x[1],x[2],_,_,_ = AA.T
    for nn in range(6):
        if (not nn==0) and (not nn==2):
            cell = cells[0,nn]
            templen = len(cell[0,0])
            trigarr = np.ones(templen)
            for k in range(templen):
                if cell[0,2][k]==1:
                    trigarr[k] = np.cos(np.asarray(cell[0,3][k]).item() * np.asarray(AA[4]).item())
                if cell[0,2][k]==2:
                    trigarr[k] = np.sin(np.asarray(cell[0,3][k]).item() * np.asarray(AA[4]).item())
            result[nn] = np.sum((cell[0,0].T * np.prod((x ** cell[0,1]),axis=1)) * trigarr)
    return result

def AApartials(x,Lpt=1,nftype='Birkhoff',method=0):
    if (not isinstance(x,np.ndarray)):
        x = np.array(x,dtype=complex)
    if nftype=='Resonant':
        return resAApartials(x,NFdata[nftype][lpts[Lpt-1]]['AApartialscell'])
    else:
        partials = NFdata[nftype][lpts[Lpt-1]]['AApartialscell']
        result = np.zeros(6)
        for n in range(3):
            cell = partials[0,n][0]
            temp = celleval(x,cell)
            result[n+3] = temp[0].real
        return result

def qptoRTB(x,V,T1,C,a,mu): #change input to a
    return V@(T1@C@x+np.array([a-mu,0,0,0,a-mu,0]))

def RTBtoqp(x,Cinv,T1inv,Vinv,a,mu): #change input to a
    return Cinv@T1inv@Vinv@x - Cinv@T1inv@np.array([a-mu,0,0,0,a-mu,0])

def numtildetransform(x,GenFuncEOMs,flag,temparr,opts):
    state0 = x
    for n in temparr:
        rv = solve_ivp(lambda t,x: sixdeval(x,GenFuncEOMs[0,n]),(0,flag),state0,t_eval=[flag],**opts).y.T
        state0 = rv[0,:]
    return state0


def AAtoRTB(AA,Lpt=1,nftype='Birkhoff',method='anl'):
    if (not isinstance(AA,np.ndarray)):
        AA = np.array(AA,dtype=complex)
    AAsize = AA.shape
    tr = False
    if len(AAsize)==1:
        if AAsize[0]==6:
            numAA = 1
        else:
            print('ERROR: Input to AAtoRTB is the wrong size. Must be (n,6), (6,n), or (6,)')
    else:
        if AAsize[1]==6:
            numAA = AAsize[0]
        else:
            tr = True
            AA = AA.T
            numAA = AAsize[1]

    data = NFdata[nftype][lpts[Lpt-1]]

    if method=='num':
        if numAA==1:
            return qptoRTB(numtildetransform(AAtotilde(AA,nftype).real,data['GenFuncEOMs'],1,temparr,optsnum),
                       V,data['T1'],data['C'],data['aLpt'],mu)
        else:
            RTB = np.zeros([numAA,6])
            for n in range(numAA):
                RTB[n,:] = qptoRTB(numtildetransform(AAtotilde(AA[n,:],nftype).real,data['GenFuncEOMs'],1,temparr,optsnum),
                        V,data['T1'],data['C'],data['aLpt'],mu)
            if tr:
                return RTB.T
            else:
                return RTB
    else:
        if numAA==1:
            return qptoRTB(sixdeval(AAtotilde(AA,nftype),data['anlqpNtoqp0']),V,data['T1'],data['C'],data['aLpt'],mu)
        else:
            RTB = np.zeros([numAA,6])
            for n in range(numAA):
                RTB[n,:] = qptoRTB(sixdeval(AAtotilde(AA[n,:],nftype),data['anlqpNtoqp0']),V,data['T1'],data['C'],data['aLpt'],mu)
            if tr:
                return RTB.T
            else:
                return RTB

def RTBtoAA(RTB,Lpt=1,nftype='Birkhoff',method='anl'):
    if (not isinstance(RTB,np.ndarray)):
        RTB = np.array(RTB)
    RTBsize = RTB.shape
    tr = False
    if len(RTBsize)==1:
        if RTBsize[0]==6:
            numRTB = 1
        else:
            print('ERROR: Input to AAtoRTB is the wrong size. Must be (n,6), (6,n), or (6,)')
    else:
        if RTBsize[1]==6:
            numRTB = RTBsize[0]
        else:
            tr = True
            RTB = RTB.T
            numRTB = RTBsize[1]

    data = NFdata[nftype][lpts[Lpt-1]]

    if method=='num':
        if numRTB==1:
            return tildetoAA(numtildetransform(RTBtoqp(RTB,data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),
                                            data['GenFuncEOMs'],-1,temparr1,optsnum),nftype)
        else:
            AA = np.zeros([numRTB,6],dtype=complex)
            for n in range(numRTB):
                AA[n,:] = tildetoAA(numtildetransform(RTBtoqp(RTB[n,:],data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),
                                            data['GenFuncEOMs'],-1,temparr1,optsnum),nftype)
            if tr:
                return AA.T
            else:
                return AA
    else:
        if numRTB==1:
            return tildetoAA(sixdeval(RTBtoqp(RTB,data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),data['anlqp0toqpN']).real,nftype)
        else:
            AA = np.zeros([numRTB,6],dtype=complex)
            for n in range(numRTB):
                AA[n,:] = tildetoAA(sixdeval(RTBtoqp(RTB[n,:],data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),data['anlqp0toqpN']).real,nftype)
            if tr:
                return AA.T
            else:
                return AA
            
def NFtoRTB(NF,Lpt=1,nftype='Birkhoff',method='anl'):
    if (not isinstance(NF,np.ndarray)):
        NF = np.array(NF,dtype=complex)
    NFsize = NF.shape
    tr = False
    if len(NFsize)==1:
        if NFsize[0]==6:
            numNF = 1
        else:
            print('ERROR: Input to NFtoRTB is the wrong size. Must be (n,6), (6,n), or (6,)')
    else:
        if NFsize[1]==6:
            numNF = NFsize[0]
        else:
            tr = True
            NF = NF.T
            numNF = NFsize[1]

    data = NFdata[nftype][lpts[Lpt-1]]

    if method=='num':
        if numNF==1:
            return qptoRTB(numtildetransform(NF,data['GenFuncEOMs'],1,temparr,optsnum),
                       V,data['T1'],data['C'],data['aLpt'],mu)
        else:
            RTB = np.zeros([numNF,6])
            for n in range(numNF):
                RTB[n,:] = qptoRTB(numtildetransform(NF[n,:],data['GenFuncEOMs'],1,temparr,optsnum),
                        V,data['T1'],data['C'],data['aLpt'],mu)
            if tr:
                return RTB.T
            else:
                return RTB
    else:
        if numNF==1:
            return qptoRTB(sixdeval(NF,data['anlqpNtoqp0']),V,data['T1'],data['C'],data['aLpt'],mu)
        else:
            RTB = np.zeros([numNF,6])
            for n in range(numNF):
                RTB[n,:] = qptoRTB(sixdeval(NF[n,:],data['anlqpNtoqp0']),V,data['T1'],data['C'],data['aLpt'],mu)
            if tr:
                return RTB.T
            else:
                return RTB
            
def RTBtoNF(RTB,Lpt=1,nftype='Birkhoff',method='anl'):
    if (not isinstance(RTB,np.ndarray)):
        RTB = np.array(RTB)
    RTBsize = RTB.shape
    tr = False
    if len(RTBsize)==1:
        if RTBsize[0]==6:
            numRTB = 1
        else:
            print('ERROR: Input to RTBtoNF is the wrong size. Must be (n,6), (6,n), or (6,)')
    else:
        if RTBsize[1]==6:
            numRTB = RTBsize[0]
        else:
            tr = True
            RTB = RTB.T
            numRTB = RTBsize[1]

    data = NFdata[nftype][lpts[Lpt-1]]

    if method=='num':
        if numRTB==1:
            return numtildetransform(RTBtoqp(RTB,data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),
                                            data['GenFuncEOMs'],-1,temparr1,optsnum)
        else:
            NF = np.zeros([numRTB,6])
            for n in range(numRTB):
                NF[n,:] = numtildetransform(RTBtoqp(RTB[n,:],data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),
                                            data['GenFuncEOMs'],-1,temparr1,optsnum)
            if tr:
                return NF.T
            else:
                return NF
    else:
        if numRTB==1:
            return sixdeval(RTBtoqp(RTB,data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),data['anlqp0toqpN'])
        else:
            NF = np.zeros([numRTB,6])
            for n in range(numRTB):
                NF[n,:] = sixdeval(RTBtoqp(RTB[n,:],data['Cinv'],data['T1inv'],Vinv,data['aLpt'],mu),data['anlqp0toqpN'])
            if tr:
                return NF.T
            else:
                return NF
            
def AAtoNF(AA,Lpt=1,nftype='Birkhoff',method=0):
    if (not isinstance(AA,np.ndarray)):
        AA = np.array(AA,dtype=complex)
    AAsize = AA.shape
    tr = False
    if len(AAsize)==1:
        if AAsize[0]==6:
            numAA = 1
        else:
            print('ERROR: Input to AAtoNF is the wrong size. Must be (n,6), (6,n), or (6,)')
    else:
        if AAsize[1]==6:
            numAA = AAsize[0]
        else:
            tr = True
            AA = AA.T
            numAA = AAsize[1]

    data = NFdata[nftype][lpts[Lpt-1]]

    if numAA==1:
        return AAtotilde(AA,nftype)
    else:
        NF = np.zeros([numAA,6])
        for n in range(numAA):
            NF[n,:] = AAtotilde(AA[n,:],nftype)
        if tr:
            return NF.T
        else:
            return NF
        
def NFtoAA(NF,Lpt=1,nftype='Birkhoff',method='anl'):
    if (not isinstance(NF,np.ndarray)):
        NF = np.array(NF)
    NFsize = NF.shape
    tr = False
    if len(NFsize)==1:
        if NFsize[0]==6:
            numNF = 1
        else:
            print('ERROR: Input to NNtoAA is the wrong size. Must be (n,6), (6,n), or (6,)')
    else:
        if NFsize[1]==6:
            numNF = NFsize[0]
        else:
            tr = True
            NF = NF.T
            numNF = NFsize[1]

    data = NFdata[nftype][lpts[Lpt-1]]

    if numNF==1:
        return tildetoAA(NF,nftype)
    else:
        AA = np.zeros([numNF,6],dtype=complex)
        for n in range(numNF):
            AA[n,:] = tildetoAA(NF[n,:],nftype)
        if tr:
            return AA.T
        else:
            return AA
        
def Heval(AA,Lpt=1,nftype='Birkhoff',method=0):
    if (not isinstance(AA,np.ndarray)):
        AA = np.array(AA,dtype=complex)
    H = NFdata[nftype][lpts[Lpt-1]]['HAA']
    Hshape = H[0,0].shape
    x = AA[:3]
    if nftype=='Resonant':
        trigarr = np.ones([Hshape[0]])
        for n in range(Hshape[0]):
            if H[0,2][n]==1:
                trigarr[n] = np.cos(H[0,3][n]*AA[4].real)
        return np.sum((H[0,0].T * np.prod((x ** H[0,1]),axis=1)) * trigarr).real
    else:
        return np.sum((H[0,0].T * np.prod((x ** H[0,1]),axis=1))).real

    
def AAprop(tspan,AA,Lpt=1,nftype='Birkhoff',method='anl',opts=optsnum):
    if (not isinstance(AA,np.ndarray)):
        AA = np.array(AA,dtype=complex)
    if nftype=='Resonant':
        AAreal = AA.real
        AAarray = solve_ivp(lambda t,x: resAApartials(x,NFdata[nftype][lpts[Lpt-1]]['AApartialscell']),(tspan[0],tspan[-1]),AAreal,t_eval=tspan,**opts).y.T
    else:
        rates = AApartials(AA,Lpt,nftype,method)
        T = len(tspan)
        onearray = np.ones(T).T
        AAarray = np.zeros([T,6],dtype=complex)
        AAarray[:,0] = AA[0]*onearray
        AAarray[:,1] = AA[1]*onearray
        AAarray[:,2] = AA[2]*onearray
        AAarray[:,3] = rates[3]*tspan.T+AA[3]*onearray
        AAarray[:,4] = rates[4]*tspan.T+AA[4]*onearray
        AAarray[:,5] = rates[5]*tspan.T+AA[5]*onearray
    return AAarray

    
def CR3BP(t, state):
    x, y, z, xdot, ydot, zdot = state

    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)

    xdd = 2 * ydot + x - (1 - mu) * (x + mu) / r1**3 - mu * (x - (1 - mu)) / r2**3
    ydd = -2 * xdot + y - (1 - mu) * y / r1**3 - mu * y / r2**3
    zdd = -(1 - mu) * z / r1**3 - mu * z / r2**3

    return np.array([xdot, ydot, zdot, xdd, ydd, zdd])

def RTBprop(tspan, state,options=opts):
    return solve_ivp(CR3BP,(tspan[0],tspan[-1]),state,t_eval=tspan,**options).y

