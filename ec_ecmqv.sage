import hashlib
import hmac 
import binascii 

#Ellitic Curves - ECMQV
#Implementation of adding and multiplying points on EllipticCurve
#Class below represents point on an EC: y^2 = x^3 + a*x + b over a finite field with 'p' elements
#Each point has its own coordinates (x,y,z)
#A point with coordinates (0,1,0) represents the point at infinity
#The method printPoint() prints point's coordinates
#The method returnCurveParameters() returns parameters of EC
class PointOnCurve:
    def __init__(self, x, y, z, a, b, p):
        F=GF(p)
        self.a=a
        self.b=b
        self.p=p
        self.x=F(x)
        self.y=F(y)
        self.z=F(z)
    def printPoint(self):
        print("(",self.x,":",self.y,":",self.z,")")
    def returnCurveParameters(self):
        return self.a, self.b, self.p
#Function addPoints(P,Q) performs point addition on an elliptic curve (EC). 
#The function takes two points P and Q as input, where each point is represented with coordinates (x, y, z) on the EC
#If the points P and Q belong to different curves (as determined by their curve parameters), a message is printed and no addition is performed
#If both points P and Q are at the point at infinity (z=0), the result is the point at infinity (0, 1, 0) on the same EC
#If one of the points (either P or Q) is at the point at infinity (z=0), the result is the other point
#If both points P and Q are not at the point at infinity (z=1), the function performs point addition as follows:
# - If P and Q have the same (x, y) coordinates, the function calculates the tangent line and finds the third point of intersection
# - If P and Q are distinct points, the function calculates the slope of the line through P and Q and finds the third point of intersection
#The result is a new point R on the same EC, and it is returned as the output of the function.
def addPoints(P, Q):
    if(P.returnCurveParameters()!=Q.returnCurveParameters()):
        print("Points from diffrent curves")
        return
    if(P.z==1 and Q.z==1):
        if(P.x==Q.x and P.y==-Q.y):
            R=PointOnCurve(0,1,0,P.a,P.b,P.p)
        else:
            if(P.x==Q.x and P.y==Q.y):
                lamb=(3*(P.x^2)+P.a)/(2*P.y)
            else:
                lamb=(Q.y-P.y)/(Q.x-P.x)
            x3=lamb^2-P.x-Q.x
            y3=lamb*(P.x-x3)-P.y
            R=PointOnCurve(x3,y3,1,P.a,P.b,P.p)
    elif(P.z==0 and Q.z==0):
        R=PointOnCurve(0,1,0,P.a,P.b,P.p)
    else:
        if(P.z==1):
            R=PointOnCurve(P.x,P.y,P.z,P.a,P.b,P.p)
        else:
            R=PointOnCurve(Q.x,Q.y,Q.z,P.a,P.b,P.p)
    return R
#Function montgomeryLadder(P, n) performs point multiplication
#It takes a point P (represented by (x, y, z) coordinates on an elliptic curve), and a scalar n as input
#It initializes two points, P0 and P1, with P0 set to the point at infinity (0, 1, 0) and P1 set to the input point P
#The binary representation of the scalar n is obtained, and its length (number of bits) is stored in 'd'
# The function then iterates through the binary representation of 'n' from the most significant bit to the least significant bit
# For each bit:
# - If the bit is 1, it performs point addition by adding 'P0' and 'P1' to it, and also doubles 'P1'
# - If the bit is 0, it performs point addition by adding 'P1' and 'P0' to it, and also doubles 'P0'
#After all iterations, the final result is in 'P0', which represents the scalar multiplication of 'P' by 'n' on the elliptic curve
#The resulting point 'P0' is returned as the output of the function
def montgomeryLadder(P, n):
    P0=PointOnCurve(0,1,0,P.a,P.b,P.p)
    P1=PointOnCurve(P.x,P.y,P.z,P.a,P.b,P.p)
    n_bit=bin(n)[2:]
    d=len(n_bit)
    for i in range (d-1,-1,-1):
        if(int(n_bit[d-1-i])==1):
            P0=addPoints(P0,P1)
            P1=addPoints(P1,P1)
        else:
            P1=addPoints(P0,P1)
            P0=addPoints(P0,P0)
    return P0


# Dash function for extracting bits from the x-coordinate of a point on an elliptic curve
# Parameters:
# - R - point on the elliptic curve
# - r - order of the elliptic curve
# Returns:
# - x - extracted bits from the x-coordinate of the poin
def dashFunction(R,r):
    f=floor(log(r,2))+1
    x=R.x
    x=int(x) % 2**(ceil(f/2))
    return x
   
# Key derivation function using PBKDF2-HMAC-SHA256
# Parameters:
# - x - input value for key derivation.
# Returns:
# - two derived keys, k1 and k2.
# Function converts x to binary representation, then encode x as bytes,
# performs key derivation using PBKDF2-HMAC-SHA256 and converts the derived key
# to hexadecimal representation. Finally splits the key into two parts, k1 and k2
def keyDerive(x):
    x=bin(x)[2:]
    a=x.encode('utf-8')
    key=hashlib.pbkdf2_hmac('sha256',a,a,1000,32)
    key=binascii.hexlify(key).decode('utf-8')
    k1=key[0:32]
    k2=key[32:]
    return k1,k2
    
# Function MAC() Compute a Message Authentication Code (MAC) using HMAC-SHA256
# Parameters:
# - nr - message number
# - P1, P2, R1, R2 - two public keys and two random public keys
# - pKey - secret key for HMAC
# Returns:
# - mac - hexadecimal representation of the computed MAC.
# Function concatenates message components then encodes the message as bytes.
# Nwxt step is computing HMAC-SHA256 using the provided secret key and
# converting the MAC to a hexadecimal string
def MAC(nr,P1,P2,R1,R2,pKey):
    mess=str(nr)+str(P1)+str(P2)+str(R1)+str(R2)
    mess=mess.encode('utf-8')
    mac = hmac.new(pKey.encode('utf-8'), mess, hashlib.sha256).digest()
    mac=mac.hex()
    return mac
    

    
# MQV function performs Elliptic Curve MQV (ECMQV) key agreement protocol
# Parameters:
# - p - ,odulus of the elliptic curve
# - a, b - coefficients of the elliptic curve equation
# - r - order of the base point G
# - Gx, Gy, Gz - coordinates of the base point G on the elliptic curve
# - h - co-factor of the elliptic curve
# Returns:
# - k2 or -1 - session key k2 if the protocol is successful, -1 otherwise.
# This function simulates the ECMQV key agreement protocol, allowing two parties (A and B) to establish
# a shared session key k2. The protocol involves private and public keys, random values, and cryptographic
# operations. If the protocol completes successfully, the derived session key k2 is returned; otherwise, -1 is returned.    
def ECMQV(p,a,b,r,Gx,Gy,Gz,h):
    
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    
    print("Elliptic Curve MQV")
    print("Domain parameters:")
    print("\tp =",p)
    print("\ta =",a)
    print("\tb =",b)
    print("\tr =",r)
    print("\th =",h)
    print("\tBase point G = ", end="")
    G.printPoint()
    
    kA=randint(2,r)
    QA=montgomeryLadder(G,kA)
    
    kB=randint(2,r)
    QB=montgomeryLadder(G,kB)
    
    print("\nA keys:")
    print("\tPrivate key - kA =", kA)
    print("\tPublic key - QA =", end="")
    QA.printPoint()
    print("\nB keys:")
    print("\tPrivate key - kB =", kB)
    print("\tPublic key - QB =", end="")
    QB.printPoint()
    
    
    dA=randint(2,r-2)
    RA=montgomeryLadder(G,dA)

    
    print("\nProtocol steps:")
    print("1. A chooses randomly dA from [2,r-2] and calculates RA=[dA]G\n\tdA =",dA,"\n\tRA =",end="")
    RA.printPoint()
    
    print("A ----- QA, RA ------> B")
    
    
    dB=randint(2,r-2)
    RB=montgomeryLadder(G,dB)
    RbDash=dashFunction(RB,r)
    sB=(dB+RbDash*kB) % r
    Z=montgomeryLadder((addPoints(RA,montgomeryLadder(QA,dashFunction(RA,r)))),(h*sB))
    x=Z.x
    x=int(x)
    k1,k2=keyDerive(x)
    tB=MAC(2,QB,QA,RB,RA,k1)

    
    print("\n2. B does:")
    print("\t1. Verifies RA")
    print("\t2. Chooses randomly dB from [2,r-2] and calculates RB=[dB]G\n\t\tdB =",dB,"\n\t\tRB =",end="")
    RB.printPoint()
    print("\t3. Calculates sB=(dB+(RB^)kB) mod r and Z=[hsB](RA+[RA^]QA), checks if Z is not point at infinity")
    print("\t\tsB =",sB)
    print("\t\tZ =",end="")
    Z.printPoint()
    if(int(Z.z)==1):
        print("\t\tZ is not point at infinity")
    else:
        print("\t\tZ is point at infinity")
        return -1
    print("\t4. Calculates (k1,k2):=KDF(xz), where xz is x coordinate of point Z")
    print("\t\tk1 =",k1)
    print("\t\tk2 =",k2)
    print("\t5. Calculates tB=MAC_k1 (2,QB,QA,RB,RA)")
    print("\t\ttB =", tB)
    
    print("\nA ----- QB, RB, tB ----> B")
    
    

    sA=(dA+dashFunction(RA,r)*kA) % r
    aZ=montgomeryLadder((addPoints(RB,(montgomeryLadder(QB,(dashFunction(RB,r)))))),(h*sA))
    ak1,ak2=keyDerive(int(aZ.x))
    t=MAC(2,QB,QA,RB,RA,ak1)
    tA=MAC(3,QA,QB,RA,RB,ak1)

    
    print("\n3. A does:")
    print("\t1. Verifies RB")
    print("\t2. Calculates sA=(dA+(RA^)kA) mod r and Z=[hsA](RB+[RB^]QB) and checks if Z is point at infinity")
    print("\t\tsA =",sA)
    print("\t\tZ =",end="")
    aZ.printPoint()
    if(int(aZ.z)==1):
        print("\t\tZ is not point at infinity")
    else:
        print("\t\tZ is point at infinity")
        return -1
    print("\t3. Calculates (k1,k2):=KDF(xz), where xz is x coordinate of point Z")
    print("\t\tk1=",ak1)
    print("\t\tk2=",ak2)
    print("\t4. Calculates t = MAC_k1 (2,QB,QA,RB,RA) and checks if t == tB")
    print("\t\tt =", t)
    if(t==tB):
        print("\t\tt equals tB")
    else:
        print("\t\tt doesn't equal tB")
        return -1
    print("\t5. Calculates tA = MAC_k1 (3,QA,QB,RA,RB)")
    print("\t\ttA =", tA)
    
    print("\nA --- tA --->B")
    
    
    bt=MAC(3,QA,QB,RA,RB,ak1)
    
    print("\n4. B calculates t = MAC_k1 (3,QA,QB,RA,RB) and checks if t equals tA")
    print("\tt =",bt)
    if(bt==tA):
        print("\tt equals tA")
    else:
        print("\tt doesn't equal tA")
        return -1
 
    print("\n5. If everything went properly k2 is session key:")
    print("\tk2 =",k2)
    return k2
    
p=39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319
a=-3
b=27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575
r=39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643
Gx=0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7
Gy=0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f 
Gz=1
G=PointOnCurve(Gx,Gy,Gz,a,b,p)        
h=1
k=ECMQV(p,a,b,r,Gx,Gy,Gz,h)