import numpy as np
import matplotlib.pyplot as plt
import random
import itertools
from itertools import combinations
from itertools import permutations


class charge:
    def __init__(self, q, pos):
        self.q=q
        self.pos=pos


"""Definting multiplication between two elements"""
def mult(p1, p2):
    """Multiplication between two elements"""
    return p1*p2


"""Return the distance between a set of charges """
def distance(p1,p2):
    """Euclidean distance between two points."""
    x1,y1 = p1
    x2,y2 = p2
    return np.hypot(x2 - x1, y2 - y1)


"""Return the force created by each set of charges """
k = 8.988 * (10**9) # In Nm^2/C^2                                                               # k is the constant in electrostatics
def F(Q, d):                                                                                                   # charges are in nC                                
    return (k*Q)/d**2


def E(q, r0, x, y):
    """Return the electric field vector E=(Ex,Ey) due to charge q at r0."""
    den = np.hypot(x-r0[0], y-r0[1])**3
    return q * (x - r0[0]) / den, q * (y - r0[1]) / den


def V_point_charge(q, r0, x, y):
    """Return the electric potential scaler due to charge q at r0."""
    d = np.hypot(((x - r0[0])*10**(-2)), ((y - r0[1])*10**(-2)))
    return k*(q*10**(-6)/d)



"""Grid of x, y points"""
nx, ny = 424, 424
axis_limit = 2
x = np.linspace(-axis_limit, axis_limit, nx)
y = np.linspace(-axis_limit, axis_limit, ny)
X, Y = np.meshgrid(x, y)


num_charges = input("What would you like to generate? Monopole (M), Dipole (D), Quadropole (Q), Octopole (O), Random (R): ")

"""Setting charges locations and values"""

#"""Generate five random charges"""
if num_charges == 'R':
    charges = []

    rq1 = np.random.uniform(-10, 10)
    rq2 = np.random.uniform(-10, 10)
    rq3 = np.random.uniform(-10, 10)
    rq4 = np.random.uniform(-10, 10)
    rq5 = np.random.uniform(-10, 10)

    rx, ry = np.random.uniform(-2, 2), np.random.uniform(-2, 2)
    rx2, ry2 = np.random.uniform(-2, 2), np.random.uniform(-2, 2)
    rx3, ry3 = np.random.uniform(-2, 2), np.random.uniform(-2, 2)
    rx4, ry4 = np.random.uniform(-2, 2), np.random.uniform(-2, 2)
    rx5, ry5 = np.random.uniform(-2, 2), np.random.uniform(-2, 2)

    charges.append((round(rq1, 3), [round(rx, 3), round(ry, 3)]))
    charges.append((round(rq2, 3), [round(rx2, 3), round(ry2, 3)]))
    charges.append((round(rq3, 3), [round(rx3, 3), round(ry3, 3)]))
    charges.append((round(rq4, 3), [round(rx4, 3), round(ry4, 3)]))
    charges.append((round(rq5, 3), [round(rx5, 3), round(ry5, 3)]))


#"""Monopole"""
elif num_charges == 'M':
    rq = round(np.random.uniform(-10, 10), 3)
    charges = [(rq, [0, 0])]

#"""DIpole"""
elif num_charges ==  'D':
    rq = round(np.random.uniform(-10, 10), 3)
    charges = [(rq, [-1, 0]), (-(rq), [1, 0])]

#"""Quadropole"""
elif num_charges ==  'Q':
    rq = round(np.random.uniform(-10, 10), 3)
    charges = [(rq, [-1, -1]), (-(rq), [1, -1]), (rq, [1, 1]), (-(rq), [-1, 1])]

#"""Octopole"""
elif num_charges ==  'O':
    rq = round(np.random.uniform(-10, 10), 3)
    nq = 8
    charges = []
    for i in range(nq):
        Aq = i%2 * 2 - 1
        q = Aq*rq
        charges.append((q, (round((2/np.sqrt(2))*np.cos(2*np.pi*i/nq),3), round((2/np.sqrt(2))*np.sin(2*np.pi*i/nq),3))))
        
else:
    print("\nNo charges generated")



"""Setting charges locations and values"""
##charges = [(2, [-1, 0]), (3, [0, 0]), (1, [0, -np.sqrt(3)])]                 #potential must be 117,000 at 0 (-1, [0, np.sqrt(3)])
#charges = [(2, [-1, 0]), (3, [0, 0]), (1, [0, 3])] 
#charges = [(-8.5, [.5 , .5]), (3, [0, np.sqrt(3)]), (4, [-1, .5]), (5, [1.5, -1.5])]                    # Formatted as (charge value, [x-value, y-value])


##################################################################################################


"""Seperating the charge values, and cooridnates from the tuples"""
coordinates = []
charge_values = []
printedCharges = []

for p in charges:
    coordinates.append(p[1])
    charge_values.append(float(p[0]))
    
    mu =  u"\u03bc"
    printedCharges.append(str(float(p[0])) + " " + mu + "C")


"""Calculating the force of each charge against all other charges using Colomb's law, then calculating the sum of the x and y components of the force on each charge"""
if len(charges) == 1:
    print("\nThere are no outside forces on this charge")

else:
    for i in range(0, len(coordinates)):
        Test_Charge = []
        Test_Charge.append(i)

        relative_distances = []
        charge_multiples = []
        thetas = []
        indicies = range(len(charges))
        for j in  (set(indicies)-set(Test_Charge)):
            x_com = coordinates[i][0] - coordinates[j][0]
            y_com =  coordinates[i][1] - coordinates[j][1] 

            relative_distances.append(np.hypot(x_com, y_com)*(10**(-2)))
            test_charge = charge_values[i]
            source_charge = charge_values[j]
            charge_multiples.append(test_charge*(10**(-6))*source_charge*(10**(-6)))

            try:
                theta =  np.degrees(np.arctan(abs(y_com/x_com)))
            except:
                theta = 0
            if x_com == 0 and y_com < 0:
                theta = 270
            elif x_com == 0 and y_com > 0:
                theta= 90
            elif x_com > 0 and y_com == 0:
                theta = 0
            elif x_com < 0 and y_com == 0:
                theta = 180
            elif x_com < 0 and y_com > 0:
                theta = 180 - theta
            elif x_com < 0 and y_com < 0:
                theta = 270 - theta
            elif x_com > 0 and y_com < 0:
                theta = 360 - theta
            else:
                theta = theta

            thetas.append(theta)
            
        Forces = []
        Fx_forces = []
        Fy_forces = []
        for Q in range(len(charge_multiples)):
            Force = (k*((charge_multiples[Q])/(relative_distances[Q]**2)))
            Forces.append(Force)
            Fx_forces.append(Force*np.cos(np.deg2rad(thetas[Q])))
            Fy_forces.append(Force*np.sin(np.deg2rad(thetas[Q])))

        Fx_total = sum(Fx_forces)
        Fy_total = sum(Fy_forces)
        Net_force = np.hypot(Fx_total, Fy_total)
        
        try:
            final_theta =  np.degrees(np.arctan(abs(Fy_total/Fx_total)))
        except:
            final_theta = 0   
        if Fx_total == 0 and Fy_total < 0:
            final_theta = 270.0
        elif Fx_total ==0 and Fy_total > 0:
            final_theta = 90.0
        elif Fx_total > 0 and Fy_total == 0:
            final_theta = 0.0
        elif Fx_total < 0 and Fy_total == 0:
            final_theta= 180.0
        elif Fx_total < 0 and Fy_total > 0:
            final_theta = 180 - final_theta
        elif Fx_total < 0 and Fy_total < 0:
            final_theta = 270 - final_theta
        elif Fx_total > 0 and Fy_total < 0:
            final_theta = 360 - final_theta
        else:
            final_theta = final_theta

        deg = u"\u00b0"
        finalString_p1 = "\nTotal Net Force on Charge" + " " + str(i+1) + " " + "is: " + " " + str("{:,}".format(round(Net_force, 1))) + " " + "N"  + " at " + str("{:,}".format(round(final_theta, 1))) + deg .ljust(20)
        finalString_p2 = ("{For charge value: " + str(printedCharges[i]) + " | " + "At coordinate point: " + str(coordinates[i]) + "}")
        print (finalString_p1 + finalString_p2)

        
################################################################################################


"""Calculate and graph the Electric Potential for each point in the grid due to the sum of each charge"""
def V_total(x, y, charges):
    V = 0
    for C in charges:
        Vp = V_point_charge(C[0], C[1], x, y)
        V  = V + Vp
    return V

Z = V_total(X, Y, charges)


################################################################################################


"""Calculate the work it takes to assemble all of the charges (The electric potential energy of the system)"""
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


"""Electric field vector, E=(Ex, Ey), as separate components"""
Ex, Ey = np.zeros((ny, nx)), np.zeros((ny, nx))
for charge in charges:
    ex, ey = E(*charge,  x=X, y=Y)
    Ex += ex
    Ey += ey


#############################################################################################

        
fig = plt.figure()
ax = fig.add_subplot(111)

"""Plot the streamlines with an appropriate colormap and arrow style"""
color = 2 * np.log(np.hypot(Ex, Ey))
ax.streamplot(x, y, Ex, Ey, color=color, linewidth=1, cmap=plt.cm.inferno,
              density=1.7, arrowstyle='->', arrowsize=1.0)


"""Add filled circles for the charges themselves"""
charge_colors = {True: '#aa0000', False: '#0000aa'}
for q, pos in charges:
    ax.add_artist(plt.Circle(pos, 0.015, color=charge_colors[q>-1]))


"""Insert the graph settings"""
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
plt.show()


