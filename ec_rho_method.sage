#Ellitic Curves
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
 
#The function exEuclides(x,N) is extended Euclidean Algorithm to find the modular inverse of 'x' modulo 'N'
# Returns:
# Tuple (u, v, d): u and v are integers such that u*x + v*N = d (the greatest common divisor)
#    If d==1, then u is the modular inverse of x modulo N
#    If d>1, then x has no modular inverse modulo N
# Note:
#    The function assumes that x should be smaller than N
#    If x is greater than or equal to N, it prints a message and returns 1
def exEuclides(x,N):
    if(x>=N):
        print("x should be smaller than N")
        return 1
    A=N
    B=x
    Ua=0
    Ub=1
    M=matrix([[A],[B]])
    Q=matrix([[0,1],[1,5]])
    U=matrix([[Ua],[Ub]])
    while True:
        q=floor(M[0,0]/M[1,0])
        Q[1,1]=-q
        M=Q*M
        U=Q*U
        if M[1,0]==0:
            break
    d=M[0,0]
    u=U[0,0]
    v=(d-x*u)/N
    return u,v,d
 
#A class representing a triplet containing three components: T, b, and g    
class Triplet:
    #Constructor method for the Triplet class
    # Parameters:
    #   T: An object representing a PointOnCurve
    #   b: Integer, one component of the triplet
    #   g: Integer, another component of the triplet
    def __init__(self, T,b,g):
        self.T=T
        self.b=b
        self.g=g
    #Method to print the components of the Triplet
    # Prints:
    #   The T component by calling the printPoint() method of the T object
    #   The b component
    #   The g component   
    def printTriplet(self):
        print("T:")
        self.T.printPoint()
        print("b:",self.b)
        print("g:",self.g)

#fFunction(Q,k) is projecting onto a scalar by taking the x-coordinate of point Q modulo k.   
# Parameters:
#   Q: An object representing a PointOnCurve
#   k: Integer, the modulus for computing the result.
# Returns:
#   x-coordinate of point Q modulo k.
def fFunction(Q,k):
    a=Integer(Q.x)
    a=a%k
    return(a)


# findingCollision(p,a,b,Gx,Gy,Gz,rG,k) is used to find a collision in a given elliptic curve
# Parameters:
# - p: Prime field
# - a, b: Curve coefficients
# - Gx, Gy, Gz: Coordinates of the base point
# - rG: Order of the base point
# - k: Number of iterations
#
# The function initializes a base point 'G' on the elliptic curve and then performs a collision search using random values.
# It generates a random value 'alpha', computes 'Q = alpha * G', and creates a list 'M' to store Triplets (directions).
# Triplets consist of points 'T', and random values 'b' and 'g'. These Triplets are generated 'k' times.
# The function then enters a collision search loop, updating points and variables until a collision is found.
# Collision search:
# - The loop iterates until a collision is found
# - In each iteratation, the loop calculates 'L1' and 'L2' values using 'fFunction' for points 'R0' and 'S0'
# - Using 'L1' and 'L2' points 'R0' and 'S0' are updated by using 'addPoints' with adequate points from list 'M' (points from triplets with 'L1' and 'L2' as indexes) 
# - Values of 'a0' and 'b0' are updated by adding modulo n to them values 'b' and 'g' from list 'M' (with 'L1' as index)
# - Values of 'c0' and 'd0' are updated by adding modulo n to them values 'b' and 'g' from list 'M' (with 'L1' as index)
# - Using 'fFunction' value 'L3' is calculated for new point 'S0'
# - Using 'L3' point 'S0' is updated by using 'addPoints' with adequate point from list 'M' (point from triplet with 'L3' as index)
# - Values of 'c0' and 'd0' are updated by adding modulo n to them values 'b' and 'g' from list 'M' (with 'L3' as index)
# - If 'R0' and 'S0' have the same x-coordinate:
#    - If 'R0' and 'S0' have the same y-coordinate:
#          - It calculates values  'A = (c0 - a0) % n' and 'B = (b0 - d0) % n'
#          - Then it checks gcd(B,n), if it doesn't eqaul one, fucntion returns to searching
#          - Else it calculates value 'foundAlpha = (A*BInv) % n' where 'BInv' is inverse of B modulo n, calculated using extended Eucildean Algorithm
#          - If calculated value equals generated value the function print information about found collision
#    - If 'R0' and 'S0' have opposite y-coordinate:
#          - It calculates values  'A = (-c0 - a0) % n' and 'B = (b0 + d0) % n'
#          - Then it checks gcd(B,n), if it doesn't eqaul one, fucntion returns to searching
#          - Else it calculates value 'foundAlpha = (A*BInv) % n' where 'BInv' is inverse of B modulo n, calculated using extended Eucildean Algorithm
#          - If calculated value equals generated value the function print information about found collision
#    - Else it returns to search loop
def findingCollision(p,a,b,Gx,Gy,Gz,rG,k): 
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    n=rG
    alpha=randint(2,n-2)
    print("Random value of alpha =",alpha)
    Q=montgomeryLadder(G,alpha)
    M=[]   
    for i in range (0, k):
        bi=randint(2,n-2)
        gi=randint(2,n-2)
        Ti=addPoints(montgomeryLadder(G,bi),montgomeryLadder(Q,gi))
        record=Triplet(Ti,bi,gi)
        M.append(record)
 
    a0=randint(2,n-2)
    b0=randint(2,n-2)
    c0=a0
    d0=b0
    R0=addPoints(montgomeryLadder(G,a0),montgomeryLadder(Q,b0))
    S0=R0
    while True:
        L1=fFunction(R0,k)
        
        Tl=M[L1].T
        bl=M[L1].b
        gl=M[L1].g
        
        R0=addPoints(R0,Tl)
        a0=(a0+bl)%n
        b0=(b0+gl)%n
        
        L2=fFunction(S0,k)
        
        Tl2=M[L2].T
        bl2=M[L2].b
        gl2=M[L2].g
        
        S0=addPoints(S0,Tl2)
        c0=(c0+bl2)%n
        d0=(d0+gl2)%n
        
        L3=fFunction(S0,k)
        
        Tl3=M[L3].T
        bl3=M[L3].b
        gl3=M[L3].g
        
        S0=addPoints(S0,Tl3)
        
        c0=(c0+bl3)%n
        d0=(d0+gl3)%n
 
        if(R0.x==S0.x):
            if(R0.y==S0.y):
                A=(c0-a0)%n
                B=(b0-d0)%n
                BInv,o,gcd=exEuclides(B,n)
                if (gcd!=1):
                    continue
                foundAlpha=(A*BInv)%n
                print("Found value of alpha = ",foundAlpha)
                if(foundAlpha==alpha):
                    print("Found alpha equals generated alpha\n")
                else:
                    print("Found alpha doesn't equal generated alpha\n")
                break;
            elif(R0.y==-S0.y):  
                A=(-c0-a0)%n
                B=(b0+d0)%n
                BInv,o,gcd=exEuclides(B,n)
                if (gcd!=1):
                    continue
                foundAlpha=(A*BInv)%n
                print("Found value of alpha = ",foundAlpha)
                if(foundAlpha==alpha):
                    print("Found alpha equals generated alpha\n")
                else:
                    print("Found alpha doesn't equal generated alpha\n")
                break;
            else: 
                continue

findingCollision(1073741789,382183198,410736703,431583365,858920426,1,1073759053,10)    
findingCollision(1099511627689,937626108435,666042130277,30009621022,215563891949,1,1099512159103,10) 