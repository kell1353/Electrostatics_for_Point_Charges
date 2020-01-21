import sys
import numpy as np
import matplotlib.pyplot as plt
import random
import itertools
from itertools import combinations
from itertools import permutations

import time
start_time = time.time()

class charge:
    def __init__(self, q, pos):
        self.q=q
        self.pos=pos

"Definting multiplication between two elements"
def mult(p1, p2):
    "Multiplication between two elements"
    return p1*p2

"Return the distance between a set of charges"
def distance(p1,p2):
    "Euclidean distance between two points."
    x1,y1 = p1
    x2,y2 = p2
    return np.hypot(x2 - x1, y2 - y1)

"Return the force created by each set of charges "
k = 8.988 * (10**9) # In Nm^2/C^2                                                               # k is the constant in electrostatics
def F(Q, d):                                                                                                                                        
    return (k*Q)/d**2

def E(q, r0, x, y):
    "Return the electric field vector E=(Ex,Ey) due to charge q at r0."
    den = np.hypot(x-r0[0], y-r0[1])**3
    return q * (x - r0[0]) / den, q * (y - r0[1]) / den

def V_point_charge(q, r0, x, y):
    "Return the electric potential scaler due to charge q at r0."
    d = np.hypot(((x - r0[0])*10**(-2)), ((y - r0[1])*10**(-2)))
    return k*(q*10**(-6)/d)


#######################################################################################
"End of Definitions"
#######################################################################################

"Common Variables"
sqrt = np.sqrt
sin, cos, tan = np.sin, np.cos, np.tan

"Create a grid of points to evaluate the magnetic field at"
nx, ny = 424, 424
axis_limit = 2
x = np.linspace(-axis_limit, axis_limit, nx)
y = np.linspace(-axis_limit, axis_limit, ny)
X, Y = np.meshgrid(x, y)


options = ['M', 'D', 'Q', 'O', 'R']
num_charges = input("What would you like to generate? Monopole (M), Dipole (D), Quadropole (Q), Octopole (O), Random (R): ")

cType = num_charges.upper()
if cType not in options:
    print("\nNo charges generated")
    sys.exit()
    
"Setting charges locations and values  (rq (uC), Distance (cm) "
########## Generate five random charges ##########
charges = []
if cType == 'R':
    totalAmt = 9

    rq, rx, ry = {}, {}, {}
    for i in range(1, totalAmt + 1):
        rq['rq'+str(i)] = np.random.uniform(-10, 10)
        rx['rx'+str(i)] = np.random.uniform(-axis_limit, axis_limit)
        ry['ry'+str(i)] = np.random.uniform(-axis_limit, axis_limit)
        
    for i in range(1, totalAmt + 1):
        charges.append((round(rq['rq'+str(i)],3), [round(rx['rx'+str(i)], 3), round(ry['ry'+str(i)], 3)]))

else:
    ##"Monopole"##
    if cType == 'M': nq = 1
    ##"DIpole"##
    elif cType ==  'D': nq = 2
    ##"Quadropole"##
    elif cType ==  'Q': nq = 4
    ##"Octopole"##
    elif cType ==  'O': nq = 8

    "Create a set of alternating charges"
    rq = round(np.random.uniform(-10, 10), 3)
    for i in range(nq):
        Aq = i%2 * 2 - 1
        q = Aq*rq
        charges.append((q, (round((2/sqrt(2))*cos(2*np.pi*i/nq),3), round((2/sqrt(2))*sin(2*np.pi*i/nq),3))))
    

"Setting charges locations and values for TESTING"
##charges = [(2, [-1, 0]), (3, [0, 0]), (1, [0, -np.sqrt(3)])]                 #potential must be 117,000 at 0 (-1, [0, np.sqrt(3)])
##charges = [(rq, [-1, -1]), (-(rq), [1, -1]), (rq, [1, 1]), (-(rq), [-1, 1])] #Quadropole
##charges = [(-8.5, [.5 , .5]), (3, [0, np.sqrt(3)]), (4, [-1, .5]), (5, [1.5, -1.5])]                    # Formatted as (charge value, [x-value, y-value])


##################################################################################################
"Seperating the charge values, and cooridnates from the tuples"
############################################################################################

coordinates, charge_values, printedCharges = [], [], []
for p in charges:
    coordinates.append(p[1])
    charge_values.append(float(p[0]))
    
    mu =  u"\u03bc"
    printedCharges.append(str(float(p[0])) + " " + mu + "C")

############################################################################################
"Calculating the force of each charge against all other charges using Colomb's law,"
"then calculating the sum of the x and y components of the force on each charge"
############################################################################################

if len(charges) == 1:
    print("\nThere are no outside forces on this charge")
else:
    for i in range(0, len(coordinates)):
        Test_Charge = []
        Test_Charge.append(i)

        relative_distances = []
        charge_multiples = []
        thetas = []
        for j in  (set(range(len(charges)))-set(Test_Charge)):
            # Calculate all of the component vectors
            x_com = coordinates[i][0] - coordinates[j][0]
            y_com =  coordinates[i][1] - coordinates[j][1] 

            # Convert uC to C and cm to m to solve for Force and insert the values into seperate lists
            relative_distances.append(np.hypot(x_com, y_com)*(10**(-2)))
            charge_multiples.append(charge_values[i]*(10**(-6))*charge_values[j]*(10**(-6)))

            # Calculate theta between the charges
            theta = np.degrees(np.arctan2(y_com, x_com))
            thetas.append(theta)
            
        Fx_forces, Fy_forces = [], []
        for Q in range(len(charge_multiples)):
            # Calculate the force between the charges and break it into component forces
            Force = F((charge_multiples[Q]), (relative_distances[Q]))
            Fx_forces.append(Force*np.cos(np.deg2rad(thetas[Q])))
            Fy_forces.append(Force*np.sin(np.deg2rad(thetas[Q])))

        # Sum the component forces and calculate the Net force on the charge
        Fx_total = sum(Fx_forces)
        Fy_total = sum(Fy_forces)
        Net_force = np.hypot(Fx_total, Fy_total)
        final_theta = np.degrees(np.arctan2(Fy_total, Fx_total))

        deg = u"\u00b0"
        finalString_p1 = "\nTotal Net Force on Charge" + " " + str(i+1) + " " + "is: " + " " + str("{:,}".format(round(Net_force, 1))) + " " + "N"  + " at " + str("{:,}".format(round(final_theta, 1))) + deg .ljust(20)
        finalString_p2 = ("{For charge value: " + str(printedCharges[i]) + " | " + "At coordinate point: " + str(coordinates[i]) + "}")
        print (finalString_p1 + finalString_p2)

################################################################################################
"Calculate and graph the Electric Potential for each point in the grid due to the sum of each charge"
############################################################################################

def V_total(x, y, charges):
    V = 0
    for C in charges:
        Vp = V_point_charge(C[0], C[1], x, y)
        V  = V + Vp
    return V

Z = V_total(X, Y, charges)

################################################################################################
"Calculate the work it takes to assemble all of the charges (The electric potential energy of the system)"
############################################################################################

if len(charges) == 1:
    print("\nThere is no electrical potential energy in this system.")

else:
    unique_distances = []
    for unique_coord in itertools.combinations(coordinates, 2):
        d = distance(unique_coord[0], unique_coord[1])*(10**(-2))
        unique_distances.append(d)
            
    unique_charge_mult = []
    for unique_charge in itertools.combinations(charge_values, 2):
        m = mult(unique_charge[0]*(10**(-6)), unique_charge[1]*(10**(-6)))
        unique_charge_mult.append(m)

    Work_Components = []
    for V in range(len(unique_charge_mult)):
            Work_Components.append((unique_charge_mult[V])/(unique_distances[V]))
            Work = k*sum(Work_Components)
            
    print("\nThe electrostatic potential energy of the system is: " + str(round(Work, 2)) + " J")

###############################################################################################
"Electric field vector, E=(Ex, Ey), as separate components"
############################################################################################

Ex, Ey = np.zeros((ny, nx)), np.zeros((ny, nx))
for charge in charges:
    ex, ey = E(*charge,  x=X, y=Y)
    Ex += ex
    Ey += ey

end_time = time.time()

############################################################################################
"Graph the field in Matplotlib"
############################################################################################

print(end_time - start_time)
start_time = time.time()

fig = plt.figure()
ax = fig.add_subplot(111)

"Plot the streamlines with an appropriate colormap and arrow style"
color = 2 * np.log(np.hypot(Ex, Ey))
ax.streamplot(x, y, Ex, Ey, color=color, linewidth=1, cmap=plt.cm.inferno,
              density=1.7, arrowstyle='->', arrowsize=1.0)


"Add filled circles for the charges themselves"
charge_colors = {True: '#aa0000', False: '#0000aa'}
for q, pos in charges:
    ax.add_artist(plt.Circle(pos, 0.015, color=charge_colors[q>-1]))


"Insert the graph settings"
Abs_charges = []
for q in charge_values:
    Abs_charges.append(abs(q))

M =  max(list(Abs_charges))

n = M*(10**7.25)
ax.set_xlabel("x" + " " + "(cm)")
ax.set_ylabel("y" + " " + "(cm)")
ax.set_xlim(-axis_limit, axis_limit)
ax.set_ylim(-axis_limit, axis_limit)
(plt.pcolor(X, Y, Z, cmap='seismic', vmin = -n, vmax = n)) #ax.add_artist #pcolor 
cbar = plt.colorbar(format = '%.0e', shrink = .75)
cbar.set_ticks(np.linspace(-n,n,10))
cbar.set_label("Electric Potential" + " " + "(Volts)")
ax.set_aspect('equal')
end_time = time.time()
print(end_time - start_time)
plt.show()

