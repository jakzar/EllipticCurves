#Ellitic Curves - arithemtic
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


#The testCurve(p, a, b, Gx, Gy, Gz, rG, t) function is used to test an elliptic curve and point multiplication using both a custom implementation and the SageMath library
#It takes several parameters: the prime field 'p', curve coefficients 'a' and 'b', base point coordinates 'Gx', 'Gy', 'Gz', the order of the base point 'rG'
#and the number of test iterations 't'.

#What the function does:
# Create a point 'G' using the given base point coordinates 'Gx', 'Gy', 'Gz', curve coefficients 'a' and 'b', and prime field 'p'
# Create a SageMath elliptic curve 'E' and a point 'S' on that curve using the provided parameters
# Print information about the curve and the base point
# Initialize a score for tracking differences between custom and SageMath implementations

# Perform 't' test iterations:
#   - Generate a random scalar 'k' between 2 and 'rG-2'
#   - Print the current iteration and 'k * G'
#   - Calculate 'T' using the custom 'montgomeryLadder' function and print the result
#   - Calculate 'L' using SageMath's point multiplication and print the result
#   - Compare 'T' and 'L' to check for differences. Update the difference score if differences are found

# Print a summary based on the difference score:
#   - If 'differenceScore' is 0, no differences were found in all iterations
#   - If 'differenceScore' is 1, differences were found in some iterations

# Print a separator for clarity.
# This function tests the custom implementation against the SageMath library to ensure the correctness of point multiplication on the given elliptic curve.
def testCurve(p,a,b,Gx,Gy,Gz,rG,t):
    #own implementation test
    G=PointOnCurve(Gx,Gy,Gz,a,b,p)
    
    #sage
    F=GF(p)
    a=F(a)
    b=F(b)
    E=EllipticCurve([a,b])
    S=E(Gx,Gy)
    
    print("Testing curve:",E)
    print("Point G:", S)
    
    differenceScore=0
    
    for i in range(t):
        k=randint(2,rG-2)
        print("\nIteration ",i+1)
        print(k,"* G")
        print("Own implementation: ")
        T=montgomeryLadder(G,k)
        T.printPoint()
        print("Sage implementation:")
        L=k*S
        print(L)
        
        if(T.x==L[0] and T.y==L[1] and T.z==L[2]):
            print("No differences")
        else:
            print("Differences")
            differenceScore=1
     
    
    print("\n\n\nSUMMARY")
    if(differenceScore==0):
        print("No differences were found in all iterations")
    else:
        print("Differences were found in some iterations")
        
    print("#############################################################################################################")

#Testing given elliptic curves    
testCurve(1073741789,382183198,410736703,431583365,858920426,1,1073759053,10)
testCurve(1099511627689,937626108435,666042130277,30009621022,215563891949,1,1099512159103,10)
testCurve(1125899906842597,514617658328474,865963734954572,559300734191994,352862582522159,1,1125899925928763,10)
testCurve(1152921504606846883,133449192748674296,309339390958398819,71033071733169680,537574381573531526,1,1152921505819822451,10)
testCurve(1180591620717411303389,984829373352706197321,1172503213559279140726,712933212623311168095,448008101342349699238,1,1180591620733222285993,10)
testCurve(1208925819614629174706111,225354284526360528563023,3764050222503696695444,296587682754061319950495,277682525629637456378157,1,1208925819614729757302087,10)
testCurve(1237940039285380274899124191,876614849831021940581906029,57055954074522725758222550,1133011678734820124222260756,913987897928984057367899527,1,1237940039285432637833556673,10)
testCurve(1267650600228229401496703205361,674423269691373715791761682647,734338331506318254542896260584,1209986394144692376486654159140,196194362409422914656559162883,1,1267650600228229857210281585449,10)

