import numpy as np
import cmath
import math

# DOING TASK 8 FIRST 
# the pre-fault voltages at all the buses:
V1_prefault = complex(1,-0.027)
V2_prefault = complex(0.978,-0.093)
V3_prefault = complex(1,0)
V4_prefault = complex(1.011,-0.174)
V5_prefault = complex(0.965,0.222)

# negative seq. impedances
zg1_neg = 0.15j
z12_neg = 0.25j
z14_neg = 0.2j
zbus4 = -0.27j
z45_neg = 0.1j
zbus5 = 0.3j
zbus2 = 0.21j
z23_neg = 0.125j
z24_neg = 0.2j
zg3_neg = 0.15j

# zero seq. impedances
zg1_0 = 0.075j
z12_0 = 1.1j
z14_0 = 0.88j
z45_0 = 0.44j
z23_0 = 0.55j
z24_0 = 0.88j
zg3_0 = 0.075j

# considering the constant impedances at bus 2, 4 and 5 in parallell 
# with the line impedances we have the following total impedances between the 
# respective buses:
z23_tot = (zbus2*z23_neg)/(zbus2 + z23_neg)
z45_tot = (zbus4*z45_neg)/(zbus4+z45_neg)

# doing the delta-star transformation for obtaining the Thevenin 
# equivalent impedances at the fault bus 5 for both the negative
# and zero seq. networks. negative first:

z1n_neg = (z12_neg*z14_neg)/(z12_neg + z14_neg + z24_neg)
z4n_neg = (z14_neg*z24_neg)/(z12_neg + z14_neg + z24_neg)
z2n_neg = (z12_neg*z24_neg)/(z12_neg + z14_neg + z24_neg)

# now the zero seq. transformation:

z1n_0 = (z12_0*z14_0)/ (z12_0+z14_0+z24_0)    
z4n_0 = (z14_0*z24_0)/ (z12_0+z14_0+z24_0)
z2n_0 = (z12_0*z24_0)/ (z12_0+z14_0+z24_0)

# and the corresponding Tvenenin impedances, first the negative:
Z_Thev_neg = zbus5+z45_tot+(((zg1_neg+z1n_neg)*(z23_tot+z2n_neg))/((zg1_neg+z1n_neg)+(z23_tot+z2n_neg)))

# and then the zero seq. Thevenin eq.
Z_Thev_0 = z45_0+z4n_0+(((zg1_0+z1n_0)*(z23_0+z2n_0))/((zg1_0+z1n_0)+(z23_0+z2n_0)))

# now to compose the negative and zero seq. Z buses, by way of Ybus neg. and zero

Ybus_neg = np.array([[(1/zg1_neg+1/z12_neg+1/z14_neg), -(1/z12_neg), 0, -(1/z14_neg), 0], 
            [-(1/z12_neg), (1/z12_neg+1/z24_neg+1/z23_tot), -(1/z23_tot), -(1/z24_neg), 0],
            [0, -(1/z23_tot), (1/z23_tot+1/zg3_neg), 0, 0],
            [-(1/z14_neg), -(1/z24_neg), 0, (1/z45_tot+1/z14_neg+z24_neg), -(1/z45_tot)],
            [0, 0, 0, -(1/z45_tot), (1/z45_tot + 1/zbus5)]])
print('Ybus_neg = ',Ybus_neg)
print('####################################################')
# and get the corresponding Zbus:
Zbus_neg = np.array(np.linalg.inv(Ybus_neg))
print('Zbus_neg =',Zbus_neg)
print('####################################################')

# now composing Ybus 0:
Ybus_0 = np.array([[(1/zg1_0+1/z14_0+1/z12_0), -(1/z12_0), 0, -(1/z14_0), 0],
          [-(1/z12_0), (1/z12_0+1/z24_0+1/z23_0), -(1/z23_0), -(1/z24_0), 0],
          [0, -(1/z23_0), (1/z23_0+1/zg3_0), 0, 0],
          [-(1/z14_0), -(1/z24_0), 0, (1/z14_0*1/z24_0+1/z45_0), -(1/z45_0)],
          [0, 0, 0, -(1/z45_0), 1/z45_0]]) 
print('Ybus_0 =',Ybus_0)
print('####################################################')

# and gte the corresponding Zbus:
Zbus_0 = np.array(np.linalg.inv(Ybus_0))
print('Zbus_0 =',Zbus_0)

# now to find the fault currents at bus 5: 
# Ik_0 = Ik_1 = Ik_2 = V5_prefault/ (Zkk_0 + Zkk_1 + Zkk_2 + 3Zf)
V5_prefault = 0.965
I5_0 = I5_1 = I5_2 = V5_prefault/ (Z_Thev_0 + 2*Z_Thev_neg)
print('I5_0 = I5_1 = I5_2 = ',I5_0)

# now to apply that [Iabc] = [A]*[I012]
a = complex(math.cos(120), math.sin(120))
a_2 =  complex(math.cos(240), math.sin(240))
print('a^2 =',a_2, 'and a =',a)
A = [[1, 1, 1], [1, a_2, a], [1, a, a_2]]
print('A =',A)

I_012 = [I5_0, I5_0, I5_0]
I_abc =  np.matmul(A,I_012)
print('I_abc =',I_abc )

# now to find the voltages at all buses:
# we use that Vi_0 = 0 - Zik_0 * Ik_0
#             Vi_1 = Vi_prefault - Zik * Ik_1
#             Vi_2 = 0 - Zik_2 * Ik_2
# where i is the voltage at bus i = 1, 2, 3, 4, 5 and k = number of faulty bus

# I am too tired and too ugly to do this smart so i just buldoze the hell
# out of it - HULK SMASH!!!

print(Zbus_0[0,4])
V1_0 = -(Zbus_0[0,4])*I5_0
#print('V1_0',V1_0) just for doublechecking
V1_1 = V1_prefault -(Zbus_neg[0,4])*I5_0
V1_2 = -(Zbus_neg[0,4] * I5_0)
V1_012 = np.array([V1_0, V1_1, V1_2])
V1_abc = np.array(np.matmul(np.linalg.inv(A), V1_012))

V2_0 = -(Zbus_0[1,4])*I5_0
V2_1 = V2_prefault -(Zbus_neg[1,4])*I5_0
V2_2 = -(Zbus_neg[1,4])*I5_0
V2_012 = np.array([V2_0, V2_1, V2_2])
V2_abc = np.array(np.matmul(np.linalg.inv(A), V2_012))

V3_0 = -(Zbus_0[2,4])*I5_0
V3_1 = V3_prefault -(Zbus_neg[2,4])*I5_0
V3_2 = -(Zbus_neg[2,4])*I5_0
V3_012 = np.array([V3_0, V3_1, V3_2])
V3_abc = np.array(np.matmul(np.linalg.inv(A), V3_012))

V4_0 = -(Zbus_0[3,4])*I5_0
V4_1 = V4_prefault -(Zbus_neg[3,4])*I5_0
V4_2 = -(Zbus_neg[3,4])*I5_0
V4_012 = np.array([V4_0, V4_1, V4_2])
V4_abc = np.array(np.matmul(np.linalg.inv(A), V4_012))

V5_0 = -(Zbus_0[4,4])*I5_0
V5_1 = V5_prefault -(Zbus_neg[4,4])*I5_0
V5_2 = -(Zbus_neg[4,4])*I5_0
V5_012 = np.array([V5_0, V5_1, V5_2])
V5_abc = np.array(np.matmul(np.linalg.inv(A), V5_012))      

# now to find the current in all phases between bus 4 and 5:
# first need to find I45_0 = (V4_0 - V5_0)/ Z45_0
                   # I45_1 = (V4_1 - V5_1)/ Z45_1
                   # I45_2 = (V4_2 - V5_2)/ Z45_2
# and then use the matrix transformation once more: [I45_abc] = [A]*[I45_012]
# NB: here we use the actual given line impedances in the data table of the project text
I45_0 = (V4_0 - V5_0)/ z45_0
I45_1 = (V4_1 - V5_1)/ z45_tot
I45_2 = (V4_2 - V5_2)/ z45_tot
I45_012 = np.array([I45_0, I45_1, I45_2])
I45_abc = np.array(np.matmul(A, I45_012))

# AND THIS IS THE END OF THE UNSYMMETRICAL SHORT CIRCUIT TASK 8

# STARTING TASK 7
# positive seq. impedances:
zg1_pos = 0.25j
z12_pos = 0.25j
z14_pos = 0.2j
z4_pos = -0.27j
z45_pos = 0.1j
z5_pos = 0.3j
z23_pos = 0.125j
z2_pos = 0.21j
z24_pos = 0.2j
zg3_pos = 0.25j

# considering the constant impedances at bus 2, 4 and 5 in parallell 
# with the line impedances we have the following total impedances between the 
# respective buses:
z23_tot_pos = (z2_pos*z23_pos)/(z2_pos + z23_pos)
z45_tot_pos = (z4_pos*z45_pos)/(z4_pos+z45_pos)

# composing the Ybus positive seq.
Ybus_pos = np.array([[(1/zg1_pos + 1/z12_pos + z14_pos), -(1/z12_pos), 0, -(1/z14_pos), 0],
                     [-(1/z12_pos), (1/z23_tot_pos + 1/z24_pos + 1/z12_pos), -(1/z23_tot_pos), -(1/z24_pos), 0],
                    [-(1/z12_pos), -(1/z23_tot_pos), (1/z23_pos + 1/zg3_pos), 0, 0], 
                    [-(1/z14_pos), (1/z24_pos), 0, (1/z45_tot_pos + 1/z14_pos + 1/z24_pos), -(1/z45_tot_pos)],
                    [0, 0, 0, -(1/z45_tot_pos), 1/z45_tot_pos ]])

# obtaining the Zbus positive seq. by inverting the Ybus
Zbus_pos = np.array(np.linalg.inv(Ybus_pos))
print('Zbus_pos = ',Zbus_pos)

# now to find the fault currents at bus 5 = I5_abc_pos = V5_prefault/ Zkk
# where Zkk = Z55 = Z_Thev_pos that must be obtained first:
z1n_pos = (z12_pos*z14_pos)/ (z12_pos+z14_pos+z24_pos)    
z4n_pos = (z14_pos*z24_pos)/ (z12_pos+z14_pos+z24_pos)
z2n_pos = (z12_pos*z24_pos)/ (z12_pos+z14_pos+z24_pos) 

Z_Thev_pos = zbus5+z45_tot_pos+(((zg1_pos+z1n_pos)*(z23_tot_pos+z2n_pos))/((zg1_pos+z1n_pos)+(z23_tot_pos+z2n_pos)))

I5_a = V5_prefault/ Z_Thev_pos
I5_b = a * I5_a
I5_c = a_2 * I5_a

# finding the voltages at all buses, same as before HULK SMASH!
# V1_a_post = V1_prefault + (-Z15/ Z55)*V1_prefault
# V2_a_post = V2_prefault + (-Z25/ Z55)*V2_prefault
# V3_a_post = V3_prefault + (-Z35/ Z55)*V3_prefault  
# V4_a_post = V4_prefault + (-Z45/ Z55)*V4_prefault  
# V5_a_post = V5_prefault - V5_prefault    

V1_a_post = V1_prefault -((Zbus_pos[0,4])/Zbus_pos[4,4])*V1_prefault
V1_b_post = a * V1_a_post
V1_c_post = a_2 * V1_a_post

V2_a_post = V2_prefault -((Zbus_pos[1,4])/Zbus_pos[4,4])*V2_prefault
V2_b_post = a * V2_a_post
V2_c_post = a_2 * V2_a_post

V3_a_post = V3_prefault -((Zbus_pos[2,4])/Zbus_pos[4,4])*V3_prefault
V3_b_post = a * V3_a_post
V3_c_post = a_2 * V3_a_post

V4_a_post = V4_prefault -((Zbus_pos[3,4])/Zbus_pos[4,4])*V4_prefault
V4_b_post = a * V4_a_post
V4_c_post = a_2 * V4_a_post

V5_a_post = V5_b_post = V5_c_post = 0

# and last but not least, finding the current flowing from bus 4 to bus 5
# first use Iij_postfault = (Vi_postfault - Vj_postfault)/ Zij - this will
# give the current in phase a and we can obtain the other b and c 
# by multiplying with a and a_2
I45_a_pos = (V4_a_post - V5_a_post)/ z45_tot_pos
I45_b_pos = a * I45_a_pos
I45_c_pos = a_2 * I45_a_pos

#DONE!



              
    
    
    