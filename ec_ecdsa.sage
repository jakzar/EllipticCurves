import hashlib
 
#Ellitic Curves - ECDSA
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

# Function genKeys() generates public and private keys for an ECDSA algorithm
# Parameters:
# - p - Modulus of the elliptic curve
# - a - Coefficient 'a' of the elliptic curve equation
# - b - Coefficient 'b' of the elliptic curve equation
# - r - Order of the base point G
# - Gx, Gy, Gz - Coordinates of the base point G on the elliptic curve
# Returns:
# - k - Private key
# - Q - Public key
def genKeys(p,a,b,r,Gx,Gy,Gz):
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    k=randint(2,r-2)
    Q=montgomeryLadder(G, k)
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


# genSign() generates a digital signature using the private key
# Parameters:
# - k - Private key
# - m - Message to be signed
# - p - Modulus of the elliptic curve
# - a - Coefficient 'a' of the elliptic curve equation
# - b - Coefficient 'b' of the elliptic curve equation
# - r - Order of the base point G
# - Gx, Gy, Gz - Coordinates of the base point G on the elliptic curve
# Returns:
# - t - Component of the digital signature
# - s - Component of the digital signature
def genSign(k,m,p,a,b,r,Gx,Gy,Gz):
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    while True:
        d=randint(2,r-2)
        T=montgomeryLadder(G, d)
        x1=T.x
        t=mod(x1, r)
        h = hashlib.new('sha384')
        h.update(m.encode())
        e=h.hexdigest()
        e=int(e, 16)
        dInv=inverse_mod(d, r)
        e=e+k*t
        s=mod(dInv*e, r)
        if(t!=0 and s!=0):
            break
    print("\n\nSignature generation:")
    print("\tInput:")
    print("\t\tm =",m)
    print("\t\tk =",k)
    print("\tAlgorithm:")
    print("\t\t1. d =",d)
    print("\t\t2.")
    print("\t\t\ta) d[G] = ", end="")
    T.printPoint()
    print("\t\t\tb) x1 =",x1)
    print("\t\t3. t =",t)
    print("\t\t4. e =",e)
    print("\t\t5. s =",s)
    print("\tOutput:")
    print("\t\tSignature(t,s) = (",t,", ", s,")")     
    return t,s

# Function verifySign() verify the digital signature using the public key
# Parameters:
# - t - Component of the digital signature
# - s - Component of the digital signature
# - Q - Public key
# - m - Message to be verified
# - p - Modulus of the elliptic curve
# - a - Coefficient 'a' of the elliptic curve equation
# - b - Coefficient 'b' of the elliptic curve equation
# - r - Order of the base point G
# - Gx, Gy, Gz - Coordinates of the base point G on the elliptic curve
def verifySign(t,s,Q,m,p,a,b,r,Gx,Gy,Gz):
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    print("Signature verification:")
    print("\tInput:")
    print("\t\tm =",m)    
    print("\t\tQ = ", end="")
    Q.printPoint()
    print("\t\t(t,s) = (",t,", ", s,")")
    
    print("\tAlgorithm:")
    print("\t\t1. Check if t, s are in the range [1, r-1]")
    if(t>r-1 or t<1):
        print("Invalid signature")
        return
    else:
        print("\t\t\ta) t is in the range [1, r-1]")
    if(s>r-1 or s<1):
        print("Invalid signature")
        return
    else:
        print("\t\t\tb) s is in the range [1, r-1]")

    h = hashlib.new('sha384')
    h.update(m.encode())
    e=h.hexdigest()
    e=int(e, 16)
    print("\t\t2. e =",e)
    s=int(s)
    w=inverse_mod(s, r)
    u1=mod(e*w, r)
    u2=mod(t*w, r)
    X=addPoints(montgomeryLadder(G,u1), montgomeryLadder(Q,u2))
    print("\t\t3. w =",w)
    print("\t\t4.")
    print("\t\t\ta) u1 =",u1)
    print("\t\t\tb) u2 =",u2)
    print("\t\t5. Point X = ", end="")
    X.printPoint()
    print("\t\t6. Check if X is not the point at infinity")
    if(X.z==0):
        print("Invalid signature")
        return
    else:
        print("\t\t\ta) X is not the point at infinity")
    x1=int(X.x)
    print("\t\t7.")
    print("\t\t\ta) x1 =",x1)
    v=mod(x1,r)
    print("\t\t\tb) v =",v)
    print("\t\t8. Check if v==t")
    if (v==t):
        print("\t\t a)",v," == ", t)
        print("\t\t b) v == t")
        print("\t\t c) Valid signature")
    else:
        print("\t\t a)",v," != ", t)
        print("\t\t b) v != t")
        print("\t\t c) Invalid signature")
 
 
p=39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319
a=-3
b=27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575
r=39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643
Gx=0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7
Gy=0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f 
Gz=1
G=PointOnCurve(Gx,Gy,Gz,a,b,p)        

k,Q=genKeys(p,a,b,r,Gx,Gy,Gz)
while True:
    ran=randint(2,r-2)
    if(ran!=k):
        break
QP=montgomeryLadder(G,ran)
m="ala ma kota"
mP="ala ma psa"
t,s=genSign(k, m,p,a,b,r,Gx,Gy,Gz)
print("\nCORRECT SIGNATURE VERIFICATION")
verifySign(t,s,Q,m,p,a,b,r,Gx,Gy,Gz)
print("\nINCORRECT SIGNATURE VERIFICATION - DIFFERENT MESSAGES")
verifySign(t,s,Q, mP,p,a,b,r,Gx,Gy,Gz)
print("\nINCORRECT SIGNATURE VERIFICATION - DIFFERENT PUBLIC KEY")
verifySign(t,s,QP, m,p,a,b,r,Gx,Gy,Gz)