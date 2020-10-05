import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# Creates the third point of one of the initial triangles which branch out from the pins
def initial_triangle(left, right, gamma):
    l_1 = np.sqrt((right[1]-left[1])**2+(right[0]-left[0])**2) # length between two known points
    l_2 = l_1 / np.tan(gamma) # length of line perpendicular to the line intersecting the two known points
    result = np.zeros(2)

    # The slope between the two known points
    slope = (right[1]-left[1]) / right[0]-left[0]
    result[0] = right[0] - l_2 * slope / np.sqrt(slope**2 + 1)
    result[1] = right[1] + l_2 / np.sqrt(slope**2 + 1)
    return result

# Given three points of a quadrilateral, return the fourth point if two of the angles are right.
def last_point(bottom, middle, top):
    result = np.zeros(2)
    xb = bottom[0]
    xm = middle[0]
    xt = top[0]
    yb = bottom[1]
    ym = middle[1]
    yt = top[1]

    # This algebra... was not fun. There tend to be removalable discontinuities which are nontrivial to factor and cancel.
    result[0] = (xt**2*(yb - ym) + (xb**2 + (yb - ym)*(yb - yt))*(ym - yt) + xm*(xt*(-yb + ym) + xb*(-ym + yt)))/(xt*(yb - ym) + xb*(ym - yt) + xm*(-yb + yt))
    result[1] = ((xm - xt)*(xb**2 + xm*xt - xb*(xm + xt) + yb*(yb - ym)) + (-xb + xm)*ym*yt + (xb - xm)*yt**2)/(xt*(-yb + ym) + xm*(yb - yt) + xb*(-ym + yt))
    return result

def generate_michell_truss(gamma, n, h):
    # These are each of the series in the truss. For n=5, there are 5, 4, 3, 2, 1 nodes in each layer (on one side plus the middle)
    first_layer = [np.array([h/(2.*np.tan(gamma/2)), 0.0])]
    left_point = (0.0, h/2)
    for i in range(1,n):
        first_layer.append(initial_triangle(left_point, first_layer[-1], gamma))
    all_layers = []
    all_layers.append(first_layer)

    for i in range(1, n):
        this_layer = []
        nodes_this_layer = n - i

        # Add node along horizontal axis
        first_nonzero_y = all_layers[i-1][1]
        on_horizon = all_layers[i-1][0]
        this_layer.append(last_point((first_nonzero_y[0], -first_nonzero_y[1]), on_horizon, first_nonzero_y))
        for j in range(1, nodes_this_layer):
            this_layer.append(last_point(this_layer[j-1], all_layers[i-1][j], all_layers[i-1][j+1]))
        all_layers.append(this_layer)
    return all_layers

def michell_truss_length(gamma, n, h):
    return generate_michell_truss(gamma, n, h)[-1][0][0]

### FOR TESTING HOW LENGTH LOOKS AS A FUNCTION OF GAMMA ###
# gamma_values = np.linspace(0.3, np.pi/2, 100)
# lengths = [michell_truss_length(gam, 1, h) for gam in gamma_values]
# plt.plot(gamma_values, lengths)
# plt.show()

### FOR FINDING GAMMA SUCH THAT L=40 FOR EACH TRUSS ###
# h=4
# for i in range(1, 6):
#     obj = lambda gam: michell_truss_length(gam, i, h) - 40
#     root = root_scalar(obj, bracket=(0.01, np.pi/2*.98), method='bisect')
#     print(root.root)


### FOR PLOTTING THE COLLECTION OF DISCRETE MICHELL TRUSSES ###
# heights = [16, 8, 4]
# for height in heights:
#     gammas = np.loadtxt('h_%i'%int(height))
#     for i, gamma in enumerate(gammas):
#         all_layers = generate_michell_truss(gamma, i+1, height)
# 
#         # Build lines for vertical struts
#         for layer in all_layers:
#             x_pts = []
#             y_pts = []
#             for point in layer:
#                 x_pts.append(point[0])
#                 y_pts.append(point[1])
#             plt.plot(x_pts, y_pts, 'b')
#             plt.plot(x_pts, -np.array(y_pts), 'b')
# 
#         # Build lines for horizontal struts
#         for layer in range(len(all_layers)):
#             x_pts = []
#             y_pts = []
#             for j in range(len(all_layers)):
#                 try:
#                     x_pts.append(all_layers[j][-layer-1][0])
#                     y_pts.append(all_layers[j][-layer-1][1])
#                 except:
#                     break
#             plt.plot(x_pts, y_pts, 'b')
#             plt.plot(x_pts, -np.array(y_pts), 'b')
# 
#         # Build lines for the initial connecting pieces
#         lay_0 = all_layers[0]
#         for point in lay_0:
#             plt.plot((0, point[0]), (height/2, point[1]), 'b')
#             plt.plot((0, point[0]), (-height/2, -point[1]), 'b')
# 
#         plt.plot(0, height/2, 'bs')
#         plt.plot(0, -height/2, 'bs')
#         plt.axis('equal')
#         plt.xlabel('x (ft)')
#         plt.xlabel('y (ft)')
#         plt.title('Michell truss w/ L=40, h=%i, gamma=%2.2f, n=%i'%(int(height), gamma, 2*(i+1)**2))
#         plt.savefig('h_%i_n%i.png'%(int(height),i+1))
#         plt.close()

### FOR CREATING INPUT FILES TO MY TRUSS ANALYSIS CODE ###
heights = [16, 8, 4]
for height in heights:
    gammas = np.loadtxt('h_%i'%int(height))
    for i_gamma, gamma in enumerate(gammas):
        all_layers = generate_michell_truss(gamma, i_gamma+1, height)

        outfile = open('h_%i_n_%i_analyze'%(height, i_gamma+1), 'w')
        outfile.write('joints\n')

        ### ADD ALL JOINTS ###

        joint_i = 0
        for layer in all_layers:
            x_pts = []
            y_pts = []
            for point in layer:
                x_pts.append(point[0])
                y_pts.append(point[1])

            # Write the points on and above the x axis
            for x, y in zip(x_pts, y_pts):
                outfile.write('%f %f\n'%(x, y))
                joint_i += 1
            # Write the points below the x axis
            for x, y in zip(x_pts[1:], y_pts[1:]):
                outfile.write('%f -%f\n'%(x, y))
                joint_i += 1

        ### ADD THE PIN AND ROLLER ###
        outfile.write('pins\n')
        outfile.write('0 %f\n'%(height/2))
        outfile.write('0 -%f\n'%(height/2))

        ### ADD ALL THE CONNECTIONS ###
        outfile.write('connections\n')

        # Start with connections from the pins to the first layer of supports
        i_pin = 0
        for i, p in enumerate(all_layers[0]):
            outfile.write('0p %ij\n'%i)
            i_pin += 1
        # Connect bottom pin to first node in first layer
        outfile.write('1p 0j\n')
        for i, p in enumerate(all_layers[0][1:]):
            outfile.write('1p %ij\n'%i_pin)
            i_pin += 1


        # create list of lists of the node indices
        layer_indices = []
        i_point = 0
        for layer in all_layers:
            this_layer_indx = []
            # Add on points above or on the x axis
            for point in layer:
                this_layer_indx.append(i_point)
                i_point += 1
            # Add on points below the x axis
            for point in layer[1:]:
                this_layer_indx.append(i_point)
                i_point += 1
            layer_indices.append(this_layer_indx)

        # Connect each layer vertically
        for layer in layer_indices:
            pts_per_side = int((len(layer)-1)/2)
            for i in range(pts_per_side):
                outfile.write('%ij %ij\n'%(layer[i], layer[i+1]))
            try:
                outfile.write('%ij %ij\n'%(layer[0], layer[pts_per_side+1]))
            except:
                pass
            for i in range(1, pts_per_side):
                outfile.write('%ij %ij\n'%(layer[pts_per_side+i], layer[pts_per_side+i+1]))

        # Connect each layer horizontally. This is more tricky.
        for i in range(len(layer_indices)-1):
            pts_per_side = int((len(layer_indices[i])-1)/2)
            pts_per_side_next = int((len(layer_indices[i+1])-1)/2)

            # outfile.write('%ij %ij\n'%(layer_indices[i][0], layer_indices[i+1][0]))
            outfile.write('%ij %ij\n'%(layer_indices[i][pts_per_side+1], layer_indices[i+1][0]))
            for k in range(len(layer_indices)-1):
                if pts_per_side_next-k>=0:
                    outfile.write('%ij %ij\n'%(layer_indices[i][pts_per_side-k], layer_indices[i+1][pts_per_side_next-k]))
            for k in range(pts_per_side+2, 2*pts_per_side+1):
                outfile.write('%ij %ij\n'%(layer_indices[i][k], layer_indices[i+1][k-2]))
                # outfile.write('%ij %ij\n'%(layer_indices[i][-k], layer_indices[i+1][-k]))
                # except:
                #     pass
            # for j in range(len(layer_indices)-1):
            #     outfile.write('%ij %ij\n'%(layer_indices[i][0], layer_indices[i+1][0]))
            # for j in range(len(layer_indices)-1):
            #     outfile.write('%ij %ij\n'%(layer_indices[i][pts_per_side], layer_indices[i+1][pts_per_side_next]))
        # TODO TODO TODO
        # for i in range(len(layer_indices)-1):
        #     try:
        #         outfile.write('%ij %ij\n'%(layer_indices[i][pts_per_side-2], layer_indices[i+1][pts_per_side-2]))
        #     except:
        #         pass

        for i_layer, layer in enumerate(layer_indices):
            pass
            # Need to manually connect first negative x node
            # to the node lying on the centerline
            # try:
            #     i_first_neg = len(layer)
            #     first_one = layer[len(all_layers[i_layer])]
            #     second_one = layer_indices[i_layer][0]
            #     outfile.write('%ij %ij\n'%(
            #         first_one,
            #         second_one
            #         ))
            # except:
            #     pass
            # for j in range(1, i_gamma+1-i_layer):
            #     try:
            #         outfile.write('%ij %ij\n'%(layer[j], layer_indices[i_layer+1][j-1]))
            #         outfile.write('%ij %ij\n'%(layer[j+int((len(layer))/2)], 
            #             layer_indices[i_layer+1][j+int((len(layer))/2)+1]))
            #     except:
            #         pass

        ### APPLY UNIT FORCE AT THE TIP ###
        outfile.write('forces\n')
        joint_i -= 1 # magic
        outfile.write('%i 0 -1\n'%joint_i)

        # Build lines for the initial connecting pieces
        # lay_0 = all_layers[0]
        # for point in lay_0:
        #     plt.plot((0, point[0]), (height/2, point[1]), 'b')
        #     plt.plot((0, point[0]), (-height/2, -point[1]), 'b')

        outfile.close()
