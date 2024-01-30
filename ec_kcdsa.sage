import hashlib
 
#Ellitic Curves - ECKCDSA
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

# hcertGen() function generates a hash-based certificate for elliptic curve domain parameters
# Parameters:
# - p, a, b, r, Gx, Gy, Gz - Elliptic curve domain parameters.
# Returns:
# - e - hexadecimal representation of the generated hash-based certificate
# This function generates a hash-based certificate for elliptic curve domain parameters using SHA-384.
# The certificate includes information about the modulus (p), coefficients (a, b), order of the base point (r),
# and coordinates of the base point (Gx, Gy, Gz).
def hcertGen(p,a,b,r,Gx,Gy,Gz):
    m="p="+str(p)+":a="+str(a)+":b="+str(b)+":r="+str(r)+":G()=("+str(Gx)+":"+str(Gy)+":"+str(Gz)+")"
    h = hashlib.new('sha384')
    h.update(m.encode())
    e=h.hexdigest()
    return e


# ECKCDSA key pair generation
# Parameters:
# - p - modulus of the elliptic curve
# - a, b - coefficients of the elliptic curve equation
# - r - order of the base point G
# - Gx, Gy, Gz - coordinates of the base point G on the elliptic curve.
# Returns:
# - k,Q - private key k and public key Q
def ECKCDSAgenKeys(p,a,b,r,Gx,Gy,Gz):
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    k=randint(2,r-2)
    kInv=inverse_mod(k,r)
    Q=montgomeryLadder(G, kInv)
    
    print("Domain parameters:")
    print("\tp =",p)
    print("\ta =",a)
    print("\tb =",b)
    print("\tr =",r)
    print("\tBase point G = ", end="")
    G.printPoint()
    print("\nGenerated keys:")
    print("\tk =", k)
    print("\tPublic key Q = ", end="")
    Q.printPoint()
    return k, Q

# ECKCDSA signature generation
# Parameters:
# - k - private key
# - m - message to be signed
# - p, a, b, r, Gx, Gy, Gz - elliptic curve domain parameters
# Returns:
# - signature (t, s)
def ECKCDSAgenSign(k,m,p,a,b,r,Gx,Gy,Gz):
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    
    while True:
        wasBigger=0
        d=randint(2,r-2)
        T=montgomeryLadder(G, d)
        x1=T.x
        x1=int(x1)
        x1=x1.to_bytes((x1.bit_length() + 7) // 8, byteorder='big')
        h = hashlib.new('sha384')
        h.update(x1)
        t=h.hexdigest()
        
        
        h = hashlib.new('sha384')
        hcert=hcertGen(p,a,b,r,Gx,Gy,Gz)
        h.update((hcert+m).encode())
        e=h.hexdigest()
        w=int(t,16) ^^ int(e,16)
        if (w>=r):
            wasBigger=1
            w=w-r
        s=mod(k*(d-w), r)
        if(s!=0):
            break
            
    print("\nSignature generation:")
    print("\t1. Randomly choosing from [2,r-2]")
    print("\t\td =",d)
    print("\t2. Calculating [d]G = (x1,y1)")
    print("\t\t[d]G =", end="")
    T.printPoint()
    print("\t3. Calculating t = H(x1)")
    print("\t\tt =",t)
    print("\t4. Calculating e = H(hcert,m)")
    print("\t\thcert =",hcert)
    print("\t\tm =",m)
    print("\t\te =",e)
    print("\t5. Calculating w = t XOR e and converitng do int")
    if(wasBigger==1):
        print("\t\tw =",w+r)
    else:
        print("\t\tw =",w)
    print("\t6. If w >= r then w = w-r")
    if(wasBigger==1):
        print("\t\tw was bigger than r")
        print("\t\tNew w = w - r =",w)
    else:
        print("\t\tw wasn't bigger than r")
        print("\t\tw =",w)
    print("\t7. Calculating s = k(d - w) mod r. If s==0 return to step 1")
    print("\t\ts =",s)
    if(s!=0):
        print("\t\ts != 0")
    else:
        print("\t\ts == 0, returning to step 1")
    print("\t8. Returning (t,s):")
    print("\t\t(",t,", ",s,")")
    return t,s
    
# ECKCDSA signature verification
# Parameters:
# - t, s - signature components
# - Q - public key
# - m - message to be verified
# - p, a, b, r, Gx, Gy, Gz - elliptic curve domain parameters
# Returns:
# - 0 if the signature is valid, -1 otherwise
def ECKCDSAverifySign(t,s,Q,m,p,a,b,r,Gx,Gy,Gz): 
    hcert=hcertGen(p,a,b,r,Gx,Gy,Gz)
    
    h = hashlib.new('sha384')
    h.update((hcert+m).encode())
    e=h.hexdigest()
    
    w=int(t,16)^^int(e,16)
    wasBigger=0
    if(w>=r):
        wasBigger=1
        w=w-r
    s=int(s)
    X=addPoints(montgomeryLadder(Q,s),montgomeryLadder(G,w))
    x1=X.x
    x1=int(x1)
    x1=x1.to_bytes((x1.bit_length() + 7) // 8, byteorder='big')
    h = hashlib.new('sha384')
    h.update(x1)
    v=h.hexdigest()

    
    print("\nSignature verification:")
    print("\t1. Checking if length of t is equal or lesser than IH bits and if s if from [r,1-r] and is an integer")
    print("\t\tt length (in bits) =", int(t,16).bit_length())
    print("\t\tIH bits = 384")
    if(int(t,16).bit_length()<=384):
        print("\t\tt length is correct")
    else:
        print("\t\tt length is incorrect")
        print("\t\tInvalid signature")
        return -1
    if(s<1 or s>r-1 or isinstance(s,int)==false):
        print("\t\ts is incorrect")
        print("\t\tInvalid signature")
        return -1
    else:
        print("\t\ts is correct")
    print("\t2. Calculating e = H(hcert,m)")
    print("\t\thcert =",hcert)
    print("\t\tm =",m)
    print("\t\te =",e)
    print("\t3. Calculating w = t XOR e and converitng do int")
    if(wasBigger==1):
        print("\t\tw =",w+r)
    else:
        print("\t\tw =",w)
    print("\t4. If w >= r then w = w-r")
    if(wasBigger==1):
        print("\t\tw was bigger than r")
        print("\t\tNew w = w - r =",w)
    else:
        print("\t\tw wasn't bigger than r")
        print("\t\tw =",w)
    print("\t5. Calculating X = [s]Q + [w]G")
    print("\t\tX =", end="")
    X.printPoint()
    print("\t6. Calculating v = H(x1)")
    print("\t\tv =",v)
    print("\t7. Checking if v == t")
    print("\t\tv =", v)
    print("\t\tt =", t)
    if(v==t):
        print("\t\tv==t - Signature is valid")
    else:
        print("\t\tSignature is invalid")
        return -1
    
    
p=39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319
a=-3
b=27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575
r=39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643
Gx=0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7
Gy=0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f 
Gz=1
G=PointOnCurve(Gx,Gy,Gz,a,b,p)        
m="ala ma kota"

k,Q=ECKCDSAgenKeys(p,a,b,r,Gx,Gy,Gz)
t,s=ECKCDSAgenSign(k,m,p,a,b,r,Gx,Gy,Gz)
ECKCDSAverifySign(t,s,Q,m,p,a,b,r,Gx,Gy,Gz)
