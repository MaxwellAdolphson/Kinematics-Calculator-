'''
Maxwell Adolphson 
Python Final Project
TA: Ranganathan Chidambaranathan
'''
#Imported to use the absolute value and the squre root finction
from math import fabs, sqrt

class PhysicsCalcKine(object):
    #Constructor, intilizes all the variables to be used in this class
    def __init__(self):
        self.fx = 0
        self.vartosolve = 0
        self.v0 = 0
        self.v = 0
        self.a = 0
        self.t = 0
        self.x = 0
        self.x0 = 0
        self.results = 0
    #Setter, to acces the function choice in main
    def setfx(self):
        return int(self.fx)
    #Briefs the user on what this part of the function can do. The two options are the ftwo functions listed below.
    def options(self):
        print ("This function will prefrom some basic Kinematic functions commonly found in introductory physics courses.")
        print ("---------------------------------------------------------------------------------------------------------")
        print("(1) V = V0 + at    or     (2) X = x0 + vt +(1/2)at^2  ")
        print ("---------------------------------------------------------------------------------------------------------")
    #Takes in the user choice of which function they want to solve.
        self.fx = raw_input("Please enter in 1 or 2 to select the function you want to solve:: ")
        print ("---------------------------------------------------------------------------------------------------------")
        self.fx= int(self.fx)

        if(self.fx ==1):
    #If the user selected the first function, they are asked which variable they would like to solve for. 
            print ("Which variable would you like to solve for?")
            print ("---------------------------------------------------------------------------------------------------------")
            self.vartosolve = raw_input("(1) V,  (2) v0,  (3) a, ( 4) t:: ")
            print ("---------------------------------------------------------------------------------------------------------")
            self.vartosolve = int(self.vartosolve)
    #The next four if statements take in the other constants, other than the one the user chose to solve for. 
            if(self.vartosolve == 1):
                print("Please enter v0, a, and t with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.v0 = raw_input()
                self.a = raw_input()
                self.t = raw_input()
                
            elif (self.vartosolve == 2):
                print("Please enter V, a, and t with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.v = raw_input()
                self.a = raw_input()
                self.t = raw_input()
                
            elif(self.vartosolve ==3):
                print("Please enter V, v0, and t with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.v = raw_input()
                self.v0 = raw_input()
                self.t = raw_input()
                
            elif (self.vartosolve == 4):
                print("Please enter V, v0, a with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.v = raw_input()
                self.v0 = raw_input()
                self.a = raw_input()
                
            


        elif(self.fx ==2):
    #If the user selected the second function, they are asked which variable they would like to solve for. 
            print ("Which variable would you like to solve for?")
            print ("There is not an option to solve for t in this function.")
            print ("---------------------------------------------------------------------------------------------------------")
            self.vartosolve = raw_input("(1) V,  (2) x0,  (3) a,  (4) X ::")
            print ("---------------------------------------------------------------------------------------------------------")
            self.vartosolve = int(self.vartosolve)
    #The next four if statements take in the other constants, other than the one the user chose to solve for.  
            if(self.vartosolve == 1):
                print("Please enter X, x0,  a, and t with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.x = raw_input()
                self.x0 = raw_input()
                self.a = raw_input()
                self.t = raw_input()  
            elif (self.vartosolve == 2):
                print("Please enter X, V, t, a with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.x = raw_input()
                self.v = raw_input()
                self.t = raw_input()
                self.a = raw_input()
            elif(self.vartosolve ==3):
                print("Please enter X, x0, v, t with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.x = raw_input()
                self.x0 = raw_input()
                self.v = raw_input()
                self.t= raw_input()
            elif (self.vartosolve == 4):
                print("Please enter x0, v, t, a with an enter inbetween")
                print ("---------------------------------------------------------------------------------------------------------")
                self.x0 = raw_input()
                self.v = raw_input()
                self.t = raw_input()
                self.a = raw_input()
    #Converts all the values taken in from the user to integers. 
        self.v =  int(self.v)
        self.v0 =  int(self.v0)
        self.a =  int(self.a)    
        self.t = int(self.t)
        self.x =  int(self.x)
        self.x0 =  int(self.x0)


    def function1(self):
    #The first function the user could have chose. This function uses the 'self.vartosolve' to access the correct version 
    #of the function, solving for the selected variable. The different 'if' statements below are different rearrangements of the 
    #Selected equation, returning self.results to be used later in the print function. 
        if(self.vartosolve == 1):
            self.results = self.v0 + (self.a*self.t)
            return self.results
            exit()
        elif(self.vartosolve == 2):
            self.results= self.v - (self.a*self.t)
            return self.results
            exit()
        elif(self.vartosolve == 3):
            self.self.results= (self.v - self.v0) / self.t
            return self.results
            exit()
        elif(self.vartosolve == 4):
            self.results= (self.v - self.v0) / self.a
            return self.results
            exit()
        else:
            return 0

    def function2(self):
    #The second function the user could have chose. This function uses the 'self.vartosolve' to access the correct version 
    #of the function, solving for the selected variable. The different 'if' statements below are different rearrangements of the 
    #Selected equation, returning self.results to be used later in the print function.
        if(self.vartosolve == 1):
            self.results = (self.x -self.x0 - (0.5)*self.a*pow(self.t,2)) / self.t
            return self.results
            exit()
        elif(self.vartosolve == 2):
            self.results= (self.x -self.v*self.t - (0.5)*self.a*pow(self.t,2))
            return self.results
            exit()
        elif(self.vartosolve == 3):
            self.results = ((2*(self.x -self.x0 - self.v*self.t) / (pow(self.t,2))))
            return self.results
            exit()
        elif(self.vartosolve == 4):
            self.results= (self.x0 +(self.v*self.t) + ((0.5)*self.a*(pow(self.t,2))))
            return self.results
            exit()
        else:
            return 0
    
    def printResults(self):
    #Print function!!!
        if(self.fx == 1):
    #If the user chose to solve the first finction, the 'self.vartosolve' is used to get into the appropriate if statment that 
    #prints "The variable you chose is "
            if(self.vartosolve == 1):
                print("The velocity is: %s"%(self.results))
            if(self.vartosolve == 2):
                print("The initial velocity is: %s"%(self.results))
            if(self.vartosolve == 3):
                print("The acceleration is: %s"%(self.results))
            if(self.vartosolve == 4):
                print("The time is: %s"%(self.results))
        elif(self.fx ==2):
    #If the user chose to solve the second finction, the 'self.vartosolve' is used to get into the appropriate if statment that 
    #prints "The variable you chose is "
            if(self.vartosolve == 1):
                print("The velocity is: %s"%(self.results))
            if(self.vartosolve == 2):
                print("The initial position is: %s"%(self.results))
            if(self.vartosolve == 3):
                print("The acceleration is: %s"%(self.results))
            if(self.vartosolve == 4):
                print("The final position is: %s"%(self.results))

class PhysicsRotCalc(PhysicsCalcKine):
    #Derived class from the PhysicsCalKine calculator
    def __inti__(self):
    #Constructor, initilizes the additional variabels to be used in the print function. 
        self.m = 0
        self.length = 0
        self.angularA = 0
        self.choice = 0
        PhysicsCalcKine.__intit__(self)

    def options(self):
    #Prompts the user about what this function does. It has the two functions from the PhysicsCalcKine class, but in rotational form. 
    #but this class also has a Torque calculator. 
        print("---------------------------------------------------------------------------------------------------------------------")
        print ("This function will prefrom some basic Rotational Kinematic functions commonly found in introductory physics courses.")
        print("---------------------------------------------------------------------------------------------------------------------")
        print("(1) w = w0 + at    or     (2) theta = theta0 + wt +(1/2)at^2  or  (3) The Torque Calculator  ")
        print("----------------------------------------------------------------------------------------------")
        self.fx = raw_input("Please enter in 1, 2 or 3 to select the function you want to solve :: ")
        print("---------------------------------------------------------------------------------")
        self.fx= int(self.fx)

        if(self.fx ==1):
    #If the user chose to solve the first function they are then asked which variable they would like to solve for. 
            print ("Which variable would you like to solve for?")
            print ("-------------------------------------------------------------------------------------------------")
            self.vartosolve = raw_input("(1) w, (2) w0, (3) a, (4) t :: ")
            self.vartosolve = int(self.vartosolve)
    #The next four 'if statements' collect the other constants the user did not chose to solve for.
            if(self.vartosolve == 1):
                print ("-------------------------------------------------------------------------------------------------")
                print("Please enter w0, a, and t with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.v0 = raw_input()
                self.a = raw_input()
                self.t = raw_input()
                
            elif (self.vartosolve == 2):
                print ("-------------------------------------------------------------------------------------------------")
                print("Please enter w, a, and t with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.v = raw_input()
                self.a = raw_input()
                self.t = raw_input()
                
            elif(self.vartosolve ==3):
                print ("-------------------------------------------------------------------------------------------------")
                print("Please enter w, w0, and t with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.v = raw_input()
                self.v0 = raw_input()
                self.t = raw_input()
                
            elif (self.vartosolve == 4):
                print ("-------------------------------------------------------------------------------------------------")
                print("Please enter w, w0, a with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.v = raw_input()
                self.v0 = raw_input()
                self.a = raw_input()
                
            


        elif(self.fx ==2):
    #If the user chose the second function they are then asked which variable they would like to  solve for. 
            print ("Which variable would you like to solve for?")
            print ("-------------------------------------------------------------------------------------------------")
            print ("There is not an option to solve for t in this function.")
            self.vartosolve = raw_input("(1) w, (2) theta0, (3) a, (4) theta ")
            print ("-------------------------------------------------------------------------------------------------")
            self.vartosolve = int(self.vartosolve)
    #The next four 'if statements' collect the other constants the user did not chose to solve for.
            if(self.vartosolve == 1):
                print("Please enter theta, theat0,  a, and t with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.x = raw_input()
                self.x0 = raw_input()
                self.a = raw_input()
                self.t = raw_input()
                exit()
            elif (self.vartosolve == 2):
                print("Please enter theta, w, t, a with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.x = raw_input()
                self.v = raw_input()
                self.t = raw_input()
                self.a = raw_input()
                exit()
            elif(self.vartosolve ==3):
                print("Please enter theta, theta0, w, t with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.x = raw_input()
                self.x0 = raw_input()
                self.v = raw_input()
                self.t= raw_input()
                exit()
            elif (self.vartosolve == 4):
                print("Please enter theta0, w, t, a with an enter inbetween")
                print ("-------------------------------------------------------------------------------------------------")
                self.x0 = raw_input()
                self.v = raw_input()
                self.t = raw_input()
                self.a = raw_input()
    #Converts all the user entered data into intergers. 
        self.v =  int(self.v)
        self.v0 =  int(self.v0)
        self.a =  int(self.a)    
        self.t = int(self.t)
        self.x =  int(self.x)
        self.x0 =  int(self.x0)

        if(self.fx ==3):
    #If the user chose to solve the torque equation. They are prompted as to what this function does and the four different shapes 
    #that are pre-programed into this function. 
            print ("This Torque calculator can calculate the torque, rotational force, of four differenct shapes.")
            print ("-------------------------------------------------------------------------------------------------")
            print ("(1) A thin rod about its center\n(2) A solid cylinder rotating about its axis.\n(3) A Solid sphere rotating about its diameter.\n(4) A hollow sphere rotating about its diameter.")
            print ("-------------------------------------------------------------------------------------------------")
            self.choice = raw_input("Please enter in your integer choice:: ")
    #Askes the user for the torque of which shape they would like to solve for. 
            print ("-------------------------------------------------------------------------------------------------")
            self.choice = int(self.choice)
    #The next four 'if statements' collect the length/radius, mass, and angular momentum of the user selected object.
            if(self.choice == 1):
                print("Please enter in the mass (kg), length (m) and angular acceleration (m/s^2) of the rod:")
                self.m = raw_input("Mass = ")
                self.length = raw_input("Radius is = ")
                self.angularA = raw_input("Angular acceleration = ")
            if(self.choice == 2):
                print("Please enter in the mass (kg), radius (m) and angular acceleration (m/s^2) of the cylinder")
                self.m = raw_input("Mass = ")
                self.length = raw_input("Radius is = ")
                self.angularA = raw_input("Angular acceleration = ")
            if(self.choice == 3):
                print("Please enter in the mass (kg), radius (m) and angular acceleration (m/s^2) of the solid sphere:")
                self.m = raw_input("Mass = ")
                self.length = raw_input("Radius is = ")
                self.angularA = raw_input("Angular acceleration = ")
            if(self.choice == 4):
                print("Please enter in the mass (kg), radius (m) and angular acceleration (m/s^2) of the hollow sphere:")
                self.m = raw_input("Mass = ")
                self.length = raw_input("Radius is = ")
                self.angularA = raw_input("Angular acceleration = ")
    def Torque(self):
    #Converts the user entered data to integers. 
        self.m = int(self.m)
        self.length = int(self.length)
        self.angularA = int(self.angularA)
    #The next four if statements are the pre-programed torque equations for the given shapes. 
        if(self.choice == 1):
            self.results = ((self.m * self.length**2) / 12) *(self.angularA)
        elif(self.choice ==2):
            self.results = ((self.m * self.length**2) / 2) *(self.angularA)
        elif(self.choice ==3):
            self.results = (2*(self.m * self.length**2) / 5) *(self.angularA)
        elif(self.choice ==4):
            self.results = (2*(self.m * self.length**2) / 3) *(self.angularA)


    def printRotResults(self):
    #This is the print function for the rotational kinematics. Although this class used the calculations of the first class
    #it needs its own print function to say angular velocity, acceleration and instead of position it is theta. 
    #Also there is a section to print the results of the torque function. 
        if(self.fx == 1):
            if(self.vartosolve == 1):
                print ("The angular velocity is %s"%(self.results))
            elif(self.vartosolve ==2):
                print ("The initial angular velocity is %s"%(self.results))
            elif (self.vartosolve == 3):
                print ("The angular acceleration is %s"%(self.results))
            elif(self.vartosolve == 4):
                print ("The time is %s"%(self.results))
        elif(self.fx ==2):
            if(self.vartosolve == 1):
                print ("The angular velocity is %s"%(self.results))
            elif(self.vartosolve == 2):
                print("The  initial angular postion is %s"%(self.results))
            elif (self.vartosolve == 3):
                print ("The angular acceleration is %s"%(self.results))
            elif(self.vartosolve ==4):
                print("The final angular position is %s"%(self.results))
        elif(self.fx == 3):
            if(self.choice == 1):
                print("The torque of the rod is: %s "%(self.results))
            elif(self.choice == 2):
                print("The torque of the solid cylinder is:  %s"%(self.results))
            elif(self.choice == 3):
                print("The torque of the solid sphere is: %s"%(self.results))
            elif(self.choice ==4):
                print("The torque of the hollow cylinder is: %s "%(self.results))
class PhysicsMathCalc(object):
    def __inti__(self):
    #Constructor, initilizes all the variables to be used in this class. 
        self.choice = 0
        self.a = 0
        self.b = 0
        self.c = 0
        self.d = 0
        self.intercept = 0
        self.upperLim = 0
        self.lowerLim = 0
        self.results = 0
        self.results2 = 0
    def options(self):
    #Pompts the user of what this part of the project does. It has three parts. 
        print("This Mathematics Calculator can prefrom three different elementary physics operations. ")
        print("---------------------------------------------------------------------------------------")
        print("(1) Newtons Method (2) Numeric Integration (3) The Quadratic Formula ")
        print("---------------------------------------------------------------------------------------")
        self.choice = raw_input("Enter in your integer choice:")
        self.choice = int(self.choice)
        if(self.choice == 1):
    #If the user selected to use Newtons method, from calculus, they are asked for the four coefficents of a cubic polynomial and a x 
    #intercept to start seraching for a zero. 
            print("This function requires the four constant coefficents of a cubic function, and the\nintercept of where you would like to search for the zero. ")
            print("i.e. Ax^3 + Bx^2 + Cx + D")
            self.a = raw_input("A= ")
            self.b = raw_input("B= ")
            self.c = raw_input("C= ")
            self.d = raw_input("D= ")
            self.intercept = raw_input("Intercept= ")
        elif(self.choice == 2):
    #If the user selected to use the numeric integration function. They are asked for the four constant coefficents of a cubic ploynomial and 
    #the upper and lower limits of integration.
            print("This program prefroms Numeric Integration of a third order polynomial using the trapzoid rule.\nThis function requires the four constant coefficents of a cubic function, and the ")
            print("the upper and lower lmits of integration")
            print("i.e. Ax^3 + Bx^2 + Cx + D")
            self.a = raw_input("A= ")
            self.b = raw_input("B= ")
            self.c = raw_input("C= ")
            self.d = raw_input("D= ")
            self.lowerLim = raw_input("Lower Limit = ")
            self.upperLim = raw_input("Upper Limit = ")
        elif(self.choice == 3):
    #If the user selected the quadratic formula, they are then asked for the three constant coefficents of a quadratic equation. 
            print("This function requires the three constant coefficients of a quadratic function\ni.e Ax^2 + Bx + C ")
            self.a = raw_input("A= ")
            self.b = raw_input("B= ")
            self.c = raw_input("C= ")
    #Converts the user entered data to integers.      
        self.a = int(self.a)
        self.b = int(self.b)
        self.c = int(self.c)
        if(self.choice == 1 or self.choice == 2):
            self.d = int(self.d)
        if(self.choice == 2):
            self.upperLim = float(self.upperLim)
            self.lowerLim = float(self.lowerLim)

    def setChoice(self):
    #Seter used in main to interact with the user. 
        return self.choice

    def NewtonsMethod(self):
    #newtons method compares the value of the derivative to the value of the original function(entered by the user), at an intercept value 
    #(entered by the user) to find the next closest value to the actual zero of the function. It uses this value in the process again. Untill 
    #the next x value and the current x value are with in a tolerance preprogrammed in to line 454 as 0.00000001
        currentx = 0.0 
        nextx=0.0 
        difference=1.0 
        dydx=0.0 
        fx=0.0
        while(difference > 0.00000001):
        #Calculating the value of the first derivative of 
        #the function at the currentx value;
            dydx = ((3*self.a*currentx**2) + 2*self.b*currentx +self.c)
        #Calculating the value of the original function at the current x value. 
            fx = (( self.a*currentx**3) + self.b*currentx**2 + self.c*currentx + self.d )
        #Finding the next closest x to the zero, in accordance with newtons method
            nextx = ((currentx)-(fx / dydx));
        #Difference tells us how close to the actual zero we are. 
            difference = fabs(nextx - currentx);
            currentx = nextx;
        self.results = currentx
    def NumericIntegration(self):
        sum_ = 0.0
        counter = 0
    #This function preforms numeric integration by means of the tapizoid meathod. This meathod requires a certain number of intervals,
    #the number of intervals equates to accuracy of this function. 
        intervals = (self.upperLim - self.lowerLim) *10000.0
    #Computes the intial x0 value and deltax value to be used in the while loop below. 
        x0 = ((self.a*self.lowerLim**3) + (self.b*self.lowerLim**2) + self.c*self.lowerLim + self.d)
        deltaX = (self.upperLim - self.lowerLim) / intervals
        
        while (counter < intervals):																																	
            counter = counter +1.0
            self.lowerLim += deltaX
            #This line sums of the area under the individual 'trapzpoids' an 'interval' number of times. 
            sum_ += 2.0 * ((self.a*self.lowerLim**3) + (self.b*self.lowerLim**2) + (self.c*self.lowerLim) + self.d )
        self.results = ((deltaX/2.0 )*(x0 + sum_)) 
    
    
    def Quadratic(self):
        discrim = 0
        notQuad = False
    #Checks to make sure the leading coefficent is not zero. 
        if(self.a == 0):
            notQuad = True
        if(notQuad):
            print("Your leading coefficient cannot be zero")

        else:
        #The quadratic formula has three different versions of output depending on the quadratic equation entered. 
        #One real, two real, or one real and one imaginary solutions. The 'if statements' below use the dicriminate to cheack for 
        #these three cases. 
            discrim = ((self.b**2) - (4*self.a *self.c))
            if(discrim < 0 ):
                discrim = sqrt(-discrim)
                discrim = -discrim
                self.results = ((-self.b)  / (2.0*self.a))
                self.results2 =  ((discrim) / (2.0*self.a))
            elif(discrim == 0):
                self.results = (-self.b) / (2*self.a)
            elif(discrim > 0):
                discrim = sqrt(discrim)
                self.results = ((-self.b + discrim) / (2.0*self.a))
                self.results2 = ((-self.b - discrim) / (2.0*self.a))


    def printResults(self):
        if(self.choice ==1):
    #Prints the results of newtons meathod. 
            print("The zero of %sx^3 + %sx^2 + %sx + %s around the point %s is located at X= %s"%(self.a, self.b, self.c, self.d, self.intercept, round(self.results,4)))
        elif(self.choice ==2):
    #Prints the results of the numeric integration function. 
            print("The approimate area under the curve %sx^3 + %sx^2 + %sx + %s is %s"%(self.a, self.b, self.c, self.d, self.results))
        elif(self.choice == 3):
    #Has three different options for the quadratic formula print funciton, listed above. 
            print("The zero of the second order poynomial %sx^2 + %sx + %s is"%(self.a, self.b, self.c))
            if(((self.b**2) - (4*self.a*self.c)) == 0):
                print(self.results)
            elif(((self.b**2) - (4*self.a*self.c)) > 0):
                print (self.results, self.results2)
            elif(((self.b**2) - (4*self.a*self.c)) < 0):
                print("%s +/- %s i"%(self.results, self.results2))

class StatisticsCalculator(object):
    def __init__(self):
    #Constructor, initilizes the variables to be used in this function. 
        self.dataArray = []
        self.sum = 0.0
        self.mean = 0.0
        self.sumSq = 0.0
    def readInData(self):
    #Reads in the numerical data from a file called "statsdata.txt" to be used in the statistics calculator. 
        results = open("statsdata.txt", "r")
        for i in range(19):
    #Assigns each line of the file to an elent of a list. 
            line = results.readline().strip()
            self.dataArray.append(line)
        
    def summation(self):
    #Returns the sum of the data
        for i in range(19):
            self.sum += int(self.dataArray[i])
        
    def Mean(self):
    #Retruns the mean value of the data
        self.mean = self.sum / 19

    def SumSquares(self):
    #Returns the sum of squares of the data. 
        for i in range(19):
            self.sumSq += (int(self.dataArray[i]))**2

    def printResults(self):
    #This print function writes the results of the calculation to a file called "statsResults.txt."
        res = open("statsResults.txt", "w")
        res.write("The sum of the data set is: %s"%(self.sum))
        res.write("\nThe mean of the data set is: %s"%(self.mean))
        res.write("\nThe sum of squares of the data set is: %s"%(self.sumSq))
        print("Data Analysis Done! Check for a file called statsResults.txt")


if __name__ == "__main__":
    #main function to operate on the classes. 
    key = True
    #While loop to continue to ask the user the below calculation choices, until the user selects the exit option. 
    while(key):
    #Prompts the user on the four parts of this project, and the exit option.
        print("\nThis program has four components. That preform some basic elementary physics Calculations, listed below")
        print("-----------------------------------------------------------------------------------------------------------------------------------")
        print("(1) [Kinematic Calculator] (2) [Rotational Kinematic Calculator] (3) [Physics Math Calculator] (4) [Satistics Calculator] (5) Exit")
        print("-----------------------------------------------------------------------------------------------------------------------------------")
        print("Please select the part of the program you would like to access for further options.")
    #The statistics calculator executes automatically when selected and immediately reads in from the file of data.
        print("If Statistics Calculator  (4) is selected the progrm will execute automatically, reading in the date from")
    #Prompts the user where to find the results of the calculation. 
        print("a file called statsdata.txt and writes the results to statsresults.txt")
        options = raw_input("Options:: ")
        print("-----------------------------------------------------------------------------------------------------------------------------------")
        options = int(options)

        
        if(options == 1):
    #Creates an object for the PhysicsCalcKine class
            k = PhysicsCalcKine()
    #Displays the options function
            k.options()
    #Uses this classes setter function to chose which function to execute. 
            if(k.setfx() == 1):
                k.function1()
            elif(k.setfx() == 2):
                k.function2()
    #Calling the print function. 
            k.printResults()

        elif(options == 2):
    #Creates an object for the PhysicsRotCalc class
            r = PhysicsRotCalc()
    #Displays the options function
            r.options()
    #Uses the setter function of the class to call the appropriate function. Or the torque calculator. 
            if(r.setfx() == 1):
                r.function1()
            elif(r.setfx == 2):
                r.function2
            elif(r.setfx() ==3):
                r.Torque()
    #Calls the print function
            r.printRotResults()

        elif(options == 3):
    #Creats and object for the PhysicsMathCalc class
            m = PhysicsMathCalc()
    #Displays the options for this class
            m.options()
    #uses the classes setetr function to all the appropriate function.
            if(m.setChoice() == 1):
                m.NewtonsMethod()
            elif(m.setChoice() ==2):
                m.NumericIntegration()  
            elif(m.setChoice() == 3):
                m.Quadratic()      
    #Calls the print function
            m.printResults()
        
        elif(options == 4):
    #Creats an object for the StatisticsCalculator. No setter function used here because this class has no options. 
            s = StatisticsCalculator()
            s.readInData()
            s.summation()
            s.Mean()
            s.SumSquares()
            s.printResults()
    #The exit choice, says goodbye and exits. 
        elif(options == 5):
            print ("GoodBye")
            key = False
