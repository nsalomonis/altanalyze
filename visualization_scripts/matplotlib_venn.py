'''
Venn diagram plotting routines.

Copyright 2012, Konstantin Tretyakov.
http://kt.era.ee/

Licensed under MIT license.

This package contains routines for plotting area-weighted two- and three-circle venn diagrams.
There are four main functions here: :code:`venn2`, :code:`venn2_circles`, :code:`venn3`, :code:`venn3_circles`.

The :code:`venn2` and :code:`venn2_circles`  accept as their only required argument a 3-element list of subset sizes:

    subsets = (Ab, aB, AB)

That is, for example, subsets[0] contains the size of the subset (A and not B), and
subsets[2] contains the size of the set (A and B), etc.

Similarly, the functions :code:`venn3` and :code:`venn3_circles` require a 7-element list:

    subsets = (Abc, aBc, ABc, abC, AbC, aBC, ABC)

The functions :code:`venn2_circles` and :code:`venn3_circles` simply draw two or three circles respectively,
such that their intersection areas correspond to the desired set intersection sizes.
Note that for a three-circle venn diagram it is not possible to achieve exact correspondence, although in
most cases the picture will still provide a decent indication.

The functions :code:`venn2` and :code:`venn3` draw diagram as a collection of separate colored patches with text labels.

The functions :code:`venn2_circles` and :code:`venn3_circles` return the list of Circle patches that may be tuned further
to your liking.

The functions :code:`venn2` and :code:`venn3` return an object of class :code:`Venn2` or :code:`Venn3` respectively,
which give access to constituent patches and text elements.

Example::

    from matplotlib import pyplot as plt
    import numpy as np
    from matplotlib_venn import venn3, venn3_circles
    plt.figure(figsize=(4,4))
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    v.get_patch_by_id('100').set_alpha(1.0)
    v.get_patch_by_id('100').set_color('white')
    v.get_label_by_id('100').set_text('Unknown')
    v.get_label_by_id('A').set_text('Set "A"')
    c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    c[0].set_lw(1.0)
    c[0].set_ls('dotted')
    plt.title("Sample Venn diagram")
    plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
                ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
'''

#### Begin __math.py

'''
Venn diagram plotting routines.
Math helper functions.

Copyright 2012, Konstantin Tretyakov.
http://kt.era.ee/

Licensed under MIT license.
'''
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

from scipy.optimize import brentq
import numpy as np

tol = 1e-10


def circle_intersection_area(r, R, d):
    '''
    Formula from: http://mathworld.wolfram.com/Circle-CircleIntersection.html
    Does not make sense for negative r, R or d

    >>> circle_intersection_area(0.0, 0.0, 0.0)
    0.0
    >>> circle_intersection_area(1.0, 1.0, 0.0)
    3.1415...
    >>> circle_intersection_area(1.0, 1.0, 1.0)
    1.2283...
    '''
    if np.abs(d) < tol:
        minR = np.min([r, R])
        return np.pi * minR**2
    if np.abs(r - 0) < tol or np.abs(R - 0) < tol:
        return 0.0
    d2, r2, R2 = float(d**2), float(r**2), float(R**2)
    arg = (d2 + r2 - R2) / 2 / d / r
    arg = np.max([np.min([arg, 1.0]), -1.0])  # Even with valid arguments, the above computation may result in things like -1.001
    A = r2 * np.arccos(arg)
    arg = (d2 + R2 - r2) / 2 / d / R
    arg = np.max([np.min([arg, 1.0]), -1.0])
    B = R2 * np.arccos(arg)
    arg = (-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R)
    arg = np.max([arg, 0])
    C = -0.5 * np.sqrt(arg)
    return A + B + C


def circle_line_intersection(center, r, a, b):
    '''
    Computes two intersection points between the circle centered at <center> and radius <r> and a line given by two points a and b.
    If no intersection exists, or if a==b, None is returned. If one intersection exists, it is repeated in the answer.

    >>> circle_line_intersection(np.array([0.0, 0.0]), 1, np.array([-1.0, 0.0]), np.array([1.0, 0.0]))
    array([[ 1.,  0.],
           [-1.,  0.]])
    >>> np.round(circle_line_intersection(np.array([1.0, 1.0]), np.sqrt(2), np.array([-1.0, 1.0]), np.array([1.0, -1.0])), 6)
    array([[ 0.,  0.],
           [ 0.,  0.]])
    '''
    s = b - a
    # Quadratic eqn coefs
    A = np.linalg.norm(s)**2
    if abs(A) < tol:
        return None
    B = 2 * np.dot(a - center, s)
    C = np.linalg.norm(a - center)**2 - r**2
    disc = B**2 - 4 * A * C
    if disc < 0.0:
        return None
    t1 = (-B + np.sqrt(disc)) / 2.0 / A
    t2 = (-B - np.sqrt(disc)) / 2.0 / A
    return np.array([a + t1 * s, a + t2 * s])


def find_distance_by_area(r, R, a, numeric_correction=0.0001):
    '''
    Solves circle_intersection_area(r, R, d) == a for d numerically (analytical solution seems to be too ugly to pursue).
    Assumes that a < pi * min(r, R)**2, will fail otherwise.

    The numeric correction parameter is used whenever the computed distance is exactly (R - r) (i.e. one circle must be inside another).
    In this case the result returned is (R-r+correction). This helps later when we position the circles and need to ensure they intersect.

    >>> find_distance_by_area(1, 1, 0, 0.0)
    2.0
    >>> round(find_distance_by_area(1, 1, 3.1415, 0.0), 4)
    0.0
    >>> d = find_distance_by_area(2, 3, 4, 0.0)
    >>> d
    3.37...
    >>> round(circle_intersection_area(2, 3, d), 10)
    4.0
    >>> find_distance_by_area(1, 2, np.pi)
    1.0001
    '''
    if r > R:
        r, R = R, r
    if np.abs(a) < tol:
        return float(r + R)
    if np.abs(min([r, R])**2 * np.pi - a) < tol:
        return np.abs(R - r + numeric_correction)
    return brentq(lambda x: circle_intersection_area(r, R, x) - a, R - r, R + r)


def circle_circle_intersection(C_a, r_a, C_b, r_b):
    '''
    Finds the coordinates of the intersection points of two circles A and B.
    Circle center coordinates C_a and C_b, should be given as tuples (or 1x2 arrays).
    Returns a 2x2 array result with result[0] being the first intersection point (to the right of the vector C_a -> C_b)
    and result[1] being the second intersection point.

    If there is a single intersection point, it is repeated in output.
    If there are no intersection points or an infinite number of those, None is returned.

    >>> circle_circle_intersection([0, 0], 1, [1, 0], 1) # Two intersection points
    array([[ 0.5      , -0.866...],
           [ 0.5      ,  0.866...]])
    >>> circle_circle_intersection([0, 0], 1, [2, 0], 1) # Single intersection point (circles touch from outside)
    array([[ 1.,  0.],
           [ 1.,  0.]])
    >>> circle_circle_intersection([0, 0], 1, [0.5, 0], 0.5) # Single intersection point (circles touch from inside)
    array([[ 1.,  0.],
           [ 1.,  0.]])
    >>> circle_circle_intersection([0, 0], 1, [0, 0], 1) is None # Infinite number of intersections (circles coincide)
    True
    >>> circle_circle_intersection([0, 0], 1, [0, 0.1], 0.8) is None # No intersections (one circle inside another)
    True
    >>> circle_circle_intersection([0, 0], 1, [2.1, 0], 1) is None # No intersections (one circle outside another)
    True
    '''
    C_a, C_b = np.array(C_a, float), np.array(C_b, float)
    v_ab = C_b - C_a
    d_ab = np.linalg.norm(v_ab)
    if np.abs(d_ab) < tol:  # No intersection points
        return None
    cos_gamma = (d_ab**2 + r_a**2 - r_b**2) / 2.0 / d_ab / r_a
    if abs(cos_gamma) > 1.0:
        return None
    sin_gamma = np.sqrt(1 - cos_gamma**2)
    u = v_ab / d_ab
    v = np.array([-u[1], u[0]])
    pt1 = C_a + r_a * cos_gamma * u - r_a * sin_gamma * v
    pt2 = C_a + r_a * cos_gamma * u + r_a * sin_gamma * v
    return np.array([pt1, pt2])


def vector_angle_in_degrees(v):
    '''
    Given a vector, returns its elevation angle in degrees (-180..180).

    >>> vector_angle_in_degrees([1, 0])
    0.0
    >>> vector_angle_in_degrees([1, 1])
    45.0
    >>> vector_angle_in_degrees([0, 1])
    90.0
    >>> vector_angle_in_degrees([-1, 1])
    135.0
    >>> vector_angle_in_degrees([-1, 0])
    180.0
    >>> vector_angle_in_degrees([-1, -1])
    -135.0
    >>> vector_angle_in_degrees([0, -1])
    -90.0
    >>> vector_angle_in_degrees([1, -1])
    -45.0
    '''
    return np.arctan2(v[1], v[0]) * 180 / np.pi


def normalize_by_center_of_mass(coords, radii):
    '''
    Given coordinates of circle centers and radii, as two arrays,
    returns new coordinates array, computed such that the center of mass of the
    three circles is (0, 0).

    >>> normalize_by_center_of_mass(np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 3.0]]), np.array([1.0, 1.0, 1.0]))
    array([[-1., -1.],
           [ 1., -1.],
           [ 0.,  2.]])
    >>> normalize_by_center_of_mass(np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]]), np.array([1.0, 1.0, np.sqrt(2.0)]))
    array([[-1., -1.],
           [ 1., -1.],
           [ 0.,  1.]])
    '''
    # Now find the center of mass.
    radii = radii**2
    sum_r = np.sum(radii)
    if sum_r < tol:
        return coords
    else:
        return coords - np.dot(radii, coords) / np.sum(radii)
        
#### End __math.py

#### Begin __venn2.py
'''
Venn diagram plotting routines.
Two-circle venn plotter.

Copyright 2012, Konstantin Tretyakov.
http://kt.era.ee/

Licensed under MIT license.
'''
import numpy as np

from matplotlib.patches import Circle
from matplotlib.colors import ColorConverter
from matplotlib.pyplot import gca

#from _math import *

#from _venn3 import make_venn3_region_patch, prepare_venn3_axes, mix_colors


def compute_venn2_areas(diagram_areas, normalize_to=1.0):
    '''
    The list of venn areas is given as 3 values, corresponding to venn diagram areas in the following order:
     (Ab, aB, AB)  (i.e. last element corresponds to the size of intersection A&B&C).
    The return value is a list of areas (A, B, AB), such that the total area is normalized
    to normalize_to. If total area was 0, returns
    (1.0, 1.0, 0.0)/2.0

    Assumes all input values are nonnegative (to be more precise, all areas are passed through and abs() function)
    >>> compute_venn2_areas((1, 1, 0))
    (0.5, 0.5, 0.0)
    >>> compute_venn2_areas((0, 0, 0))
    (0.5, 0.5, 0.0)
    >>> compute_venn2_areas((1, 1, 1), normalize_to=3)
    (2.0, 2.0, 1.0)
    >>> compute_venn2_areas((1, 2, 3), normalize_to=6)
    (4.0, 5.0, 3.0)
    '''
    # Normalize input values to sum to 1
    areas = np.array(np.abs(diagram_areas), float)
    total_area = np.sum(areas)
    if np.abs(total_area) < tol:
        return (0.5, 0.5, 0.0)
    else:
        areas = areas / total_area * normalize_to
        return (areas[0] + areas[2], areas[1] + areas[2], areas[2])


def solve_venn2_circles(venn_areas):
    '''
    Given the list of "venn areas" (as output from compute_venn2_areas, i.e. [A, B, AB]),
    finds the positions and radii of the two circles.
    The return value is a tuple (coords, radii), where coords is a 2x2 array of coordinates and
    radii is a 2x1 array of circle radii.

    Assumes the input values to be nonnegative and not all zero.
    In particular, the first two values must be positive.

    >>> c, r = solve_venn2_circles((1, 1, 0))
    >>> np.round(r, 3)
    array([ 0.564,  0.564])
    >>> c, r = solve_venn2_circles(compute_venn2_areas((1, 2, 3)))
    >>> np.round(r, 3)
    array([ 0.461,  0.515])
    '''
    (A_a, A_b, A_ab) = map(float, venn_areas)
    r_a, r_b = np.sqrt(A_a / np.pi), np.sqrt(A_b / np.pi)
    radii = np.array([r_a, r_b])
    if A_ab > tol:
        # Nonzero intersection
        coords = np.zeros((2, 2))
        coords[1][0] = find_distance_by_area(radii[0], radii[1], A_ab)
    else:
        # Zero intersection
        coords = np.zeros((2, 2))
        coords[1][0] = radii[0] + radii[1] + np.mean(radii) * 1.1
    coords = normalize_by_center_of_mass(coords, radii)
    return (coords, radii)


def compute_venn2_regions(centers, radii):
    '''
    See compute_venn3_regions for explanations.
    >>> centers, radii = solve_venn2_circles((1, 1, 0.5))
    >>> regions = compute_venn2_regions(centers, radii)
    '''
    intersection = circle_circle_intersection(centers[0], radii[0], centers[1], radii[1])
    if intersection is None:
        # Two circular regions
        regions = [("CIRCLE", (centers[a], radii[a], True), centers[a]) for a in [0, 1]] + [None]
    else:
        # Three curved regions
        regions = []
        for (a, b) in [(0, 1), (1, 0)]:
            # Make region a&not b:  [(AB, A-), (BA, B+)]
            points = np.array([intersection[a], intersection[b]])
            arcs = [(centers[a], radii[a], False), (centers[b], radii[b], True)]
            if centers[a][0] < centers[b][0]:
                # We are to the left
                label_pos_x = (centers[a][0] - radii[a] + centers[b][0] - radii[b]) / 2.0
            else:
                # We are to the right
                label_pos_x = (centers[a][0] + radii[a] + centers[b][0] + radii[b]) / 2.0
            label_pos = np.array([label_pos_x, centers[a][1]])
            regions.append((points, arcs, label_pos))

        # Make region a&b: [(AB, A+), (BA, B+)]
        (a, b) = (0, 1)
        points = np.array([intersection[a], intersection[b]])
        arcs = [(centers[a], radii[a], True), (centers[b], radii[b], True)]
        label_pos_x = (centers[a][0] + radii[a] + centers[b][0] - radii[b]) / 2.0
        label_pos = np.array([label_pos_x, centers[a][1]])
        regions.append((points, arcs, label_pos))
    return regions


def compute_venn2_colors(set_colors):
    '''
    Given two base colors, computes combinations of colors corresponding to all regions of the venn diagram.
    returns a list of 3 elements, providing colors for regions (10, 01, 11).

    >>> compute_venn2_colors(('r', 'g'))
    (array([ 1.,  0.,  0.]), array([ 0. ,  0.5,  0. ]), array([ 0.7 ,  0.35,  0.  ]))
    '''
    ccv = ColorConverter()
    base_colors = [np.array(ccv.to_rgb(c)) for c in set_colors]
    return (base_colors[0], base_colors[1], mix_colors(base_colors[0], base_colors[1]))


def venn2_circles(subsets, normalize_to=1.0, alpha=1.0, color='black', linestyle='solid', linewidth=2.0, **kwargs):
    '''
    Plots only the two circles for the corresponding Venn diagram.
    Useful for debugging or enhancing the basic venn diagram.
    parameters sets and normalize_to are the same as in venn2()
    kwargs are passed as-is to matplotlib.patches.Circle.
    returns a list of three Circle patches.

    >>> c = venn2_circles((1, 2, 3))
    >>> c = venn2_circles({'10': 1, '01': 2, '11': 3}) # Same effect
    '''
    complete_subets = subsets
    subsets_abrev = []
    for i in subsets:
        subsets_abrev.append(len(i))
    subsets = tuple(subsets_abrev)
    
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['10', '01', '11']]
    areas = compute_venn2_areas(subsets, normalize_to)
    centers, radii = solve_venn2_circles(areas)
    ax = gca()
    prepare_venn3_axes(ax, centers, radii)
    result = []
    for (c, r) in zip(centers, radii):
        circle = Circle(c, r, alpha=alpha, edgecolor=color, facecolor='none', linestyle=linestyle, linewidth=linewidth, **kwargs)
        ax.add_patch(circle)
        result.append(circle)
    return result


class Venn2:
    '''
    A container for a set of patches and patch labels and set labels, which make up the rendered venn diagram.
    '''
    id2idx = {'10': 0, '01': 1, '11': 2, 'A': 0, 'B': 1}

    def __init__(self, patches, subset_labels, set_labels):
        self.patches = patches
        self.subset_labels = subset_labels
        self.set_labels = set_labels

    def get_patch_by_id(self, id):
        '''Returns a patch by a "region id". A region id is a string '10', '01' or '11'.'''
        return self.patches[self.id2idx[id]]

    def get_label_by_id(self, id):
        '''
        Returns a subset label by a "region id". A region id is a string '10', '01' or '11'.
        Alternatively, if the string 'A' or 'B' is given, the label of the
        corresponding set is returned (or None).'''
        if len(id) == 1:
            return self.set_labels[self.id2idx[id]] if self.set_labels is not None else None
        else:
            return self.subset_labels[self.id2idx[id]]


def venn2(subsets, set_labels=('A', 'B'), set_colors=('r', 'g'), alpha=0.4, normalize_to=1.0):
    '''Plots a 2-set area-weighted Venn diagram.
    The subsets parameter is either a dict or a list.
     - If it is a dict, it must map regions to their sizes.
       The regions are identified via two-letter binary codes ('10', '01', and '11'), hence a valid set could look like:
       {'01': 10, '01': 20, '11': 40}. Unmentioned codes are considered to map to 0.
     - If it is a list, it must have 3 elements, denoting the sizes of the regions in the following order:
       (10, 10, 11)

    Set labels parameter is a list of two strings - set labels. Set it to None to disable set labels.
    The set_colors parameter should be a list of two elements, specifying the "base colors" of the two circles.
    The color of circle intersection will be computed based on those.

    The normalize_to parameter specifies the total (on-axes) area of the circles to be drawn. Sometimes tuning it (together
    with the overall fiture size) may be useful to fit the text labels better.
    The return value is a Venn2 object, that keeps references to the Text and Patch objects used on the plot.

    >>> from matplotlib_venn import *
    >>> v = venn2(subsets={'10': 1, '01': 1, '11': 1}, set_labels = ('A', 'B'))
    >>> c = venn2_circles(subsets=(1, 1, 1), linestyle='dashed')
    >>> v.get_patch_by_id('10').set_alpha(1.0)
    >>> v.get_patch_by_id('10').set_color('white')
    >>> v.get_label_by_id('10').set_text('Unknown')
    >>> v.get_label_by_id('A').set_text('Set A')
    '''
    
    complete_subets = subsets
    subsets_abrev = []
    for i in subsets:
        subsets_abrev.append(len(i))
    subsets = subsets_abrev
    
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['10', '01', '11']]
    areas = compute_venn2_areas(subsets, normalize_to)
    centers, radii = solve_venn2_circles(areas)
    if (areas[0] < tol or areas[1] < tol):
        raise Exception("Both circles in the diagram must have positive areas.")
    centers, radii = solve_venn2_circles(areas)
    regions = compute_venn2_regions(centers, radii)
    colors = compute_venn2_colors(set_colors)

    ax = gca()
    prepare_venn3_axes(ax, centers, radii)
    # Create and add patches and text
    patches = [make_venn3_region_patch(r) for r in regions]
    for (p, c) in zip(patches, colors):
        if p is not None:
            p.set_facecolor(c)
            p.set_edgecolor('none')
            p.set_alpha(alpha)
            ax.add_patch(p)
    texts = [ax.text(r[2][0], r[2][1], str(s), va='center', ha='center', size = 17) if r is not None else None for (r, s) in zip(regions, subsets)]

    # Position labels
    if set_labels is not None:
        padding = np.mean([r * 0.1 for r in radii])
        label_positions = [centers[0] + np.array([0.0, - radii[0] - padding]),
                           centers[1] + np.array([0.0, - radii[1] - padding])]
        labels = [ax.text(pos[0], pos[1], txt, size=20, ha='right', va='top') for (pos, txt) in zip(label_positions, set_labels)]
        labels[1].set_ha('left')
    else:
        labels = None
    return Venn2(patches, texts, labels)

#### End __venn2.py

#### Begin __venn3.py
'''
Venn diagram plotting routines.
Three-circle venn plotter.

Copyright 2012, Konstantin Tretyakov.
http://kt.era.ee/

Licensed under MIT license.
'''
import numpy as np
import warnings

from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path
from matplotlib.colors import ColorConverter
from matplotlib.pyplot import gca

#from _math import *

def compute_venn3_areas(diagram_areas, normalize_to=1.0):
    '''
    The list of venn areas is given as 7 values, corresponding to venn diagram areas in the following order:
     (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    (i.e. last element corresponds to the size of intersection A&B&C).
    The return value is a list of areas (A_a, A_b, A_c, A_ab, A_bc, A_ac, A_abc),
    such that the total area of all circles is normalized to normalize_to. If total area was 0, returns
    (1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)/3.0

    Assumes all input values are nonnegative (to be more precise, all areas are passed through and abs() function)
    >>> compute_venn3_areas((1, 1, 0, 1, 0, 0, 0))
    (0.33..., 0.33..., 0.33..., 0.0, 0.0, 0.0, 0.0)
    >>> compute_venn3_areas((0, 0, 0, 0, 0, 0, 0))
    (0.33..., 0.33..., 0.33..., 0.0, 0.0, 0.0, 0.0)
    >>> compute_venn3_areas((1, 1, 1, 1, 1, 1, 1), normalize_to=7)
    (4.0, 4.0, 4.0, 2.0, 2.0, 2.0, 1.0)
    >>> compute_venn3_areas((1, 2, 3, 4, 5, 6, 7), normalize_to=56/2)
    (16.0, 18.0, 22.0, 10.0, 13.0, 12.0, 7.0)
    '''
    # Normalize input values to sum to 1
    areas = np.array(np.abs(diagram_areas), float)
    total_area = np.sum(areas)
    if np.abs(total_area) < tol:
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0)
    else:
        areas = areas / total_area * normalize_to
        A_a = areas[0] + areas[2] + areas[4] + areas[6]
        A_b = areas[1] + areas[2] + areas[5] + areas[6]
        A_c = areas[3] + areas[4] + areas[5] + areas[6]

        # Areas of the three intersections (ab, ac, bc)
        A_ab, A_ac, A_bc = areas[2] + areas[6], areas[4] + areas[6], areas[5] + areas[6]

        return (A_a, A_b, A_c, A_ab, A_bc, A_ac, areas[6])


def solve_venn3_circles(venn_areas):
    '''
    Given the list of "venn areas" (as output from compute_venn3_areas, i.e. [A, B, C, AB, BC, AC, ABC]),
    finds the positions and radii of the three circles.
    The return value is a tuple (coords, radii), where coords is a 3x2 array of coordinates and
    radii is a 3x1 array of circle radii.

    Assumes the input values to be nonnegative and not all zero.
    In particular, the first three values must all be positive.

    The overall match is only approximate (to be precise, what is matched are the areas of the circles and the
    three pairwise intersections).

    >>> c, r = solve_venn3_circles((1, 1, 1, 0, 0, 0, 0))
    >>> np.round(r, 3)
    array([ 0.564,  0.564,  0.564])
    >>> c, r = solve_venn3_circles(compute_venn3_areas((1, 2, 40, 30, 4, 40, 4)))
    >>> np.round(r, 3)
    array([ 0.359,  0.476,  0.453])
    '''
    (A_a, A_b, A_c, A_ab, A_bc, A_ac, A_abc) = map(float, venn_areas)
    r_a, r_b, r_c = np.sqrt(A_a / np.pi), np.sqrt(A_b / np.pi), np.sqrt(A_c / np.pi)
    intersection_areas = [A_ab, A_bc, A_ac]
    radii = np.array([r_a, r_b, r_c])

    # Hypothetical distances between circle centers that assure
    # that their pairwise intersection areas match the requirements.
    dists = [find_distance_by_area(radii[i], radii[j], intersection_areas[i]) for (i, j) in [(0, 1), (1, 2), (2, 0)]]

    # How many intersections have nonzero area?
    num_nonzero = sum(np.array([A_ab, A_bc, A_ac]) > tol)

    # Handle four separate cases:
    #    1. All pairwise areas nonzero
    #    2. Two pairwise areas nonzero
    #    3. One pairwise area nonzero
    #    4. All pairwise areas zero.

    if num_nonzero == 3:
        # The "generic" case, simply use dists to position circles at the vertices of a triangle.
        # Before we need to ensure that resulting circles can be at all positioned on a triangle,
        # use an ad-hoc fix.
        for i in range(3):
            i, j, k = (i, (i + 1) % 3, (i + 2) % 3)
            if dists[i] > dists[j] + dists[k]:
                dists[i] = 0.8 * (dists[j] + dists[k])
                warnings.warn("Bad circle positioning")
        coords = position_venn3_circles_generic(radii, dists)
    elif num_nonzero == 2:
        # One pair of circles is not intersecting.
        # In this case we can position all three circles in a line
        # The two circles that have no intersection will be on either sides.
        for i in range(3):
            if intersection_areas[i] < tol:
                (left, right, middle) = (i, (i + 1) % 3, (i + 2) % 3)
                coords = np.zeros((3, 2))
                coords[middle][0] = dists[middle]
                coords[right][0] = dists[middle] + dists[right]
                # We want to avoid the situation where left & right still intersect
                if coords[left][0] + radii[left] > coords[right][0] - radii[right]:
                    mid = (coords[left][0] + radii[left] + coords[right][0] - radii[right]) / 2.0
                    coords[left][0] = mid - radii[left] - 1e-5
                    coords[right][0] = mid + radii[right] + 1e-5
                break
    elif num_nonzero == 1:
        # Only one pair of circles is intersecting, and one circle is independent.
        # Position all on a line first two intersecting, then the free one.
        for i in range(3):
            if intersection_areas[i] > tol:
                (left, right, side) = (i, (i + 1) % 3, (i + 2) % 3)
                coords = np.zeros((3, 2))
                coords[right][0] = dists[left]
                coords[side][0] = dists[left] + radii[right] + radii[side] * 1.1  # Pad by 10%
                break
    else:
        # All circles are non-touching. Put them all in a sequence
        coords = np.zeros((3, 2))
        coords[1][0] = radii[0] + radii[1] * 1.1
        coords[2][0] = radii[0] + radii[1] * 1.1 + radii[1] + radii[2] * 1.1

    coords = normalize_by_center_of_mass(coords, radii)
    return (coords, radii)


def position_venn3_circles_generic(radii, dists):
    '''
    Given radii = (r_a, r_b, r_c) and distances between the circles = (d_ab, d_bc, d_ac),
    finds the coordinates of the centers for the three circles so that they form a proper triangle.
    The current positioning method puts the center of A and B on a horizontal line y==0,
    and C just below.

    Returns a 3x2 array with circle center coordinates in rows.

    >>> position_venn3_circles_generic((1, 1, 1), (0, 0, 0))
    array([[ 0.,  0.],
           [ 0.,  0.],
           [ 0., -0.]])
    >>> position_venn3_circles_generic((1, 1, 1), (2, 2, 2))
    array([[ 0.        ,  0.        ],
           [ 2.        ,  0.        ],
           [ 1.        , -1.73205081]])
    '''
    (d_ab, d_bc, d_ac) = dists
    (r_a, r_b, r_c) = radii
    coords = np.array([[0, 0], [d_ab, 0], [0, 0]], float)
    C_x = (d_ac**2 - d_bc**2 + d_ab**2) / 2.0 / d_ab if np.abs(d_ab) > tol else 0.0
    C_y = -np.sqrt(d_ac**2 - C_x**2)
    coords[2, :] = C_x, C_y
    return coords


def compute_venn3_regions(centers, radii):
    '''
    Given the 3x2 matrix with circle center coordinates, and a 3-element list (or array) with circle radii [as returned from solve_venn3_circles],
    returns the 7 regions, comprising the venn diagram.
    Each region is given as [array([pt_1, pt_2, pt_3]), (arc_1, arc_2, arc_3), label_pos] where each pt_i gives the coordinates of a point,
    and each arc_i is in turn a triple (circle_center, circle_radius, direction), and label_pos is the recommended center point for
    positioning region label.

    The region is the poly-curve constructed by moving from pt_1 to pt_2 along arc_1, then to pt_3 along arc_2 and back to pt_1 along arc_3.
    Arc direction==True denotes positive (CCW) direction.

     There is also a special case, where the region is given as
    ["CIRCLE", (center, radius, True), label_pos], which corresponds to a completely circular region.

    Regions are returned in order (Abc, aBc, ABc, abC, AbC, aBC, ABC)

    >>> centers, radii = solve_venn3_circles((1, 1, 1, 1, 1, 1, 1))
    >>> regions = compute_venn3_regions(centers, radii)
    '''
    # First compute all pairwise circle intersections
    intersections = [circle_circle_intersection(centers[i], radii[i], centers[j], radii[j]) for (i, j) in [(0, 1), (1, 2), (2, 0)]]
    regions = []
    # Regions [Abc, aBc, abC]
    for i in range(3):
        (a, b, c) = (i, (i + 1) % 3, (i + 2) % 3)
        if intersections[a] is not None and intersections[c] is not None:
            # Current circle intersects both of the other circles.
            if intersections[b] is not None:
                # .. and the two other circles intersect, this is either the "normal" situation
                #    or it can also be a case of bad placement
                if np.linalg.norm(intersections[b][0] - centers[a]) < radii[a]:
                    # In the "normal" situation we use the scheme [(BA, B+), (BC, C+), (AC, A-)]
                    points = np.array([intersections[a][1], intersections[b][0], intersections[c][1]])
                    arcs = [(centers[b], radii[b], True), (centers[c], radii[c], True), (centers[a], radii[a], False)]

                    # Ad-hoc label positioning
                    pt_a = intersections[b][0]
                    pt_b = intersections[b][1]
                    pt_c = circle_line_intersection(centers[a], radii[a], pt_a, pt_b)
                    if pt_c is None:
                        label_pos = circle_circle_intersection(centers[b], radii[b] + 0.1 * radii[a], centers[c], radii[c] + 0.1 * radii[c])[0]
                    else:
                        label_pos = 0.5 * (pt_c[1] + pt_a)
                else:
                    # This is the "bad" situation (basically one disc covers two touching disks)
                    # We use the scheme [(BA, B+), (AB, A-)] if (AC is inside B) and
                    #                   [(CA, C+), (AC, A-)] otherwise
                    if np.linalg.norm(intersections[c][0] - centers[b]) < radii[b]:
                        points = np.array([intersections[a][1], intersections[a][0]])
                        arcs = [(centers[b], radii[b], True), (centers[a], radii[a], False)]
                    else:
                        points = np.array([intersections[c][0], intersections[c][1]])
                        arcs = [(centers[c], radii[c], True), (centers[a], radii[a], False)]
                    label_pos = centers[a]
            else:
                # .. and the two other circles do not intersect. This means we are in the "middle" of a OoO placement.
                # The patch is then a [(AB, B-), (BA, A+), (AC, C-), (CA, A+)]
                points = np.array([intersections[a][0], intersections[a][1], intersections[c][1], intersections[c][0]])
                arcs = [(centers[b], radii[b], False), (centers[a], radii[a], True), (centers[c], radii[c], False), (centers[a], radii[a], True)]
                # Label will be between the b and c circles
                leftc, rightc = (b, c) if centers[b][0] < centers[c][0] else (c, b)
                label_x = ((centers[leftc][0] + radii[leftc]) + (centers[rightc][0] - radii[rightc])) / 2.0
                label_y = centers[a][1] + radii[a] / 2.0
                label_pos = np.array([label_x, label_y])
        elif intersections[a] is None and intersections[c] is None:
            # Current circle is completely separate from others
            points = "CIRCLE"
            arcs = (centers[a], radii[a], True)
            label_pos = centers[a]
        else:
            # Current circle intersects one of the other circles
            other_circle = b if intersections[a] is not None else c
            other_circle_intersection = a if intersections[a] is not None else c
            i1, i2 = (0, 1) if intersections[a] is not None else (1, 0)
            # The patch is a [(AX, A-), (XA, X+)]
            points = np.array([intersections[other_circle_intersection][i1], intersections[other_circle_intersection][i2]])
            arcs = [(centers[a], radii[a], False), (centers[other_circle], radii[other_circle], True)]
            if centers[a][0] < centers[other_circle][0]:
                # We are to the left
                label_pos_x = (centers[a][0] - radii[a] + centers[other_circle][0] - radii[other_circle]) / 2.0
            else:
                # We are to the right
                label_pos_x = (centers[a][0] + radii[a] + centers[other_circle][0] + radii[other_circle]) / 2.0
            label_pos = np.array([label_pos_x, centers[a][1]])
        regions.append((points, arcs, label_pos))

    (a, b, c) = (0, 1, 2)

    # Regions [aBC, AbC, ABc]
    for i in range(3):
        (a, b, c) = (i, (i + 1) % 3, (i + 2) % 3)

        if intersections[b] is None:  # No region there
            regions.append(None)
            continue

        has_middle_region = np.linalg.norm(intersections[b][0] - centers[a]) < radii[a]

        if has_middle_region:
            # This is the "normal" situation (i.e. all three circles have a common area)
            # We then use the scheme [(CB, C+), (CA, A-), (AB, B+)]
            points = np.array([intersections[b][1], intersections[c][0], intersections[a][0]])
            arcs = [(centers[c], radii[c], True), (centers[a], radii[a], False), (centers[b], radii[b], True)]
            # Ad-hoc label positioning
            pt_a = intersections[b][1]
            dir_to_a = pt_a - centers[a]
            dir_to_a = dir_to_a / np.linalg.norm(dir_to_a)
            pt_b = centers[a] + dir_to_a * radii[a]
            label_pos = 0.5 * (pt_a + pt_b)
        else:
            # This is the situation, where there is no common area
            # Then the corresponding area is made by scheme [(CB, C+), (BC, B+), None]
            points = np.array([intersections[b][1], intersections[b][0]])
            arcs = [(centers[c], radii[c], True), (centers[b], radii[b], True)]
            label_pos = 0.5 * (intersections[b][1] + intersections[b][0])

        regions.append((points, arcs, label_pos))

    # Central region made by scheme [(BC, B+), (AB, A+), (CA, C+)]
    (a, b, c) = (0, 1, 2)
    if intersections[a] is None or intersections[b] is None or intersections[c] is None:
        # No middle region
        regions.append(None)
    else:
        points = np.array([intersections[b][0], intersections[a][0], intersections[c][0]])
        label_pos = np.mean(points, 0)  # Middle of the central region
        arcs = [(centers[b], radii[b], True), (centers[a], radii[a], True), (centers[c], radii[c], True)]
        has_middle_region = np.linalg.norm(intersections[b][0] - centers[a]) < radii[a]
        if has_middle_region:
            regions.append((points, arcs, label_pos))
        else:
            regions.append(([], [], label_pos))

    #      (Abc,        aBc,        ABc,        abC,        AbC,        aBC,        ABC)
    return (regions[0], regions[1], regions[5], regions[2], regions[4], regions[3], regions[6])


def make_venn3_region_patch(region):
    '''
    Given a venn3 region (as returned from compute_venn3_regions) produces a Patch object,
    depicting the region as a curve.

    >>> centers, radii = solve_venn3_circles((1, 1, 1, 1, 1, 1, 1))
    >>> regions = compute_venn3_regions(centers, radii)
    >>> patches = [make_venn3_region_patch(r) for r in regions]
    '''
    if region is None or len(region[0]) == 0:
        return None
    if region[0] == "CIRCLE":
        return Circle(region[1][0], region[1][1])
    pts, arcs, label_pos = region
    path = [pts[0]]
    for i in range(len(pts)):
        j = (i + 1) % len(pts)
        (center, radius, direction) = arcs[i]
        fromangle = vector_angle_in_degrees(pts[i] - center)
        toangle = vector_angle_in_degrees(pts[j] - center)
        if direction:
            vertices = Path.arc(fromangle, toangle).vertices
        else:
            vertices = Path.arc(toangle, fromangle).vertices
            vertices = vertices[np.arange(len(vertices) - 1, -1, -1)]
        vertices = vertices * radius + center
        path = path + list(vertices[1:])
    codes = [1] + [4] * (len(path) - 1)
    return PathPatch(Path(path, codes))

def mix_colors(col1, col2, col3=None):
    '''
    Mixes two colors to compute a "mixed" color (for purposes of computing
    colors of the intersection regions based on the colors of the sets.
    Note that we do not simply compute averages of given colors as those seem
    too dark for some default configurations. Thus, we lighten the combination up a bit.
    
    Inputs are (up to) three RGB triples of floats 0.0-1.0 given as numpy arrays.
    
    >>> mix_colors(np.array([1.0, 0., 0.]), np.array([1.0, 0., 0.])) # doctest: +NORMALIZE_WHITESPACE
    array([ 1.,  0.,  0.])
    >>> mix_colors(np.array([1.0, 1., 0.]), np.array([1.0, 0.9, 0.]), np.array([1.0, 0.8, 0.1])) # doctest: +NORMALIZE_WHITESPACE
    array([ 1. ,  1. , 0.04])    
    '''
    if col3 is None:
        mix_color = 0.7 * (col1 + col2)
    else:
        mix_color = 0.4 * (col1 + col2 + col3)
    mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)    
    return mix_color

def compute_venn3_colors(set_colors):
    '''
    Given three base colors, computes combinations of colors corresponding to all regions of the venn diagram.
    returns a list of 7 elements, providing colors for regions (100, 010, 110, 001, 101, 011, 111).

    >>> compute_venn3_colors(['r', 'g', 'b'])
    (array([ 1.,  0.,  0.]),..., array([ 0.4,  0.2,  0.4]))
    '''
    ccv = ColorConverter()
    base_colors = [np.array(ccv.to_rgb(c)) for c in set_colors]
    return (base_colors[0], base_colors[1], mix_colors(base_colors[0], base_colors[1]), base_colors[2],
            mix_colors(base_colors[0], base_colors[2]), mix_colors(base_colors[1], base_colors[2]), mix_colors(base_colors[0], base_colors[1], base_colors[2]))


def prepare_venn3_axes(ax, centers, radii):
    '''
    Sets properties of the axis object to suit venn plotting. I.e. hides ticks, makes proper xlim/ylim.
    '''
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    min_x = min([centers[i][0] - radii[i] for i in range(len(radii))])
    max_x = max([centers[i][0] + radii[i] for i in range(len(radii))])
    min_y = min([centers[i][1] - radii[i] for i in range(len(radii))])
    max_y = max([centers[i][1] + radii[i] for i in range(len(radii))])
    ax.set_xlim([min_x - 0.1, max_x + 0.1])
    ax.set_ylim([min_y - 0.1, max_y + 0.1])
    ax.set_axis_off()


def venn3_circles(subsets, normalize_to=1.0, alpha=1.0, color='black', linestyle='solid', linewidth=2.0, **kwargs):
    '''
    Plots only the three circles for the corresponding Venn diagram.
    Useful for debugging or enhancing the basic venn diagram.
    parameters sets and normalize_to are the same as in venn3()
    kwargs are passed as-is to matplotlib.patches.Circle.
    returns a list of three Circle patches.

    >>> plot = venn3_circles({'001': 10, '100': 20, '010': 21, '110': 13, '011': 14})
    '''
    complete_subets = subsets
    subsets_abrev = []
    for i in subsets:
        subsets_abrev.append(len(i))
    subsets = tuple(subsets_abrev)
    
    # Prepare parameters
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['100', '010', '110', '001', '101', '011', '111']]
    areas = compute_venn3_areas(subsets, normalize_to)
    centers, radii = solve_venn3_circles(areas)
    ax = gca()
    prepare_venn3_axes(ax, centers, radii)
    result = []
    for (c, r) in zip(centers, radii):
        circle = Circle(c, r, alpha=alpha, edgecolor=color, facecolor='none', linestyle=linestyle, linewidth=linewidth, **kwargs)
        ax.add_patch(circle)
        result.append(circle)
    return result


class Venn3:
    '''
    A container for a set of patches and patch labels and set labels, which make up the rendered venn diagram.
    '''
    id2idx = {'100': 0, '010': 1, '110': 2, '001': 3, '101': 4, '011': 5, '111': 6, 'A': 0, 'B': 1, 'C': 2}

    def __init__(self, patches, subset_labels, set_labels):
        self.patches = patches
        self.subset_labels = subset_labels
        self.set_labels = set_labels

    def get_patch_by_id(self, id):
        '''Returns a patch by a "region id". A region id is a string like 001, 011, 010, etc.'''
        return self.patches[self.id2idx[id]]

    def get_label_by_id(self, id):
        '''
        Returns a subset label by a "region id". A region id is a string like 001, 011, 010, etc.
        Alternatively, if you provide either of 'A', 'B' or 'C', you will obtain the label of the
        corresponding set (or None).'''
        if len(id) == 1:
            return self.set_labels[self.id2idx[id]] if self.set_labels is not None else None
        else:
            return self.subset_labels[self.id2idx[id]]


def venn3(subsets, set_labels=('A', 'B', 'C'), set_colors=('r', 'g', 'b'), alpha=0.4, normalize_to=1.0):
    global coordinates
    coordinates={}
    
    '''Plots a 3-set area-weighted Venn diagram.
    The subsets parameter is either a dict or a list.
     - If it is a dict, it must map regions to their sizes.
       The regions are identified via three-letter binary codes ('100', '010', etc), hence a valid set could look like:
       {'001': 10, '010': 20, '110':30, ...}. Unmentioned codes are considered to map to 0.
     - If it is a list, it must have 7 elements, denoting the sizes of the regions in the following order:
       (100, 010, 110, 001, 101, 011, 111).

    Set labels parameter is a list of three strings - set labels. Set it to None to disable set labels.
    The set_colors parameter should be a list of three elements, specifying the "base colors" of the three circles.
    The colors of circle intersections will be computed based on those.

    The normalize_to parameter specifies the total (on-axes) area of the circles to be drawn. Sometimes tuning it (together
    with the overall fiture size) may be useful to fit the text labels better.
    The return value is a Venn3 object, that keeps references to the Text and Patch objects used on the plot.

    >>> from matplotlib_venn import *
    >>> v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    >>> c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    >>> v.get_patch_by_id('100').set_alpha(1.0)
    >>> v.get_patch_by_id('100').set_color('white')
    >>> v.get_label_by_id('100').set_text('Unknown')
    >>> v.get_label_by_id('C').set_text('Set C')
    '''
    # Prepare parameters
    complete_subets = subsets
    subsets_abrev = []
    for i in subsets:
        subsets_abrev.append(len(i))
    subsets = tuple(subsets_abrev)
    
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['100', '010', '110', '001', '101', '011', '111']]

    areas = compute_venn3_areas(subsets, normalize_to)
    if (areas[0] < tol or areas[1] < tol or areas[2] < tol):
        raise Exception("All three circles in the diagram must have positive areas. Use venn2 or just a circle to draw diagrams with two or one circle.")
    centers, radii = solve_venn3_circles(areas)
    regions = compute_venn3_regions(centers, radii)
    colors = compute_venn3_colors(set_colors)

    ax = gca()
    prepare_venn3_axes(ax, centers, radii)
    # Create and add patches and text
    patches = [make_venn3_region_patch(r) for r in regions]
    for (p, c) in zip(patches, colors):
        if p is not None:
            p.set_facecolor(c)
            p.set_edgecolor('none')
            p.set_alpha(alpha)
            ax.add_patch(p)
    subset_labels = [ax.text(r[2][0], r[2][1], str(s), va='center', ha='center', size=17) if r is not None else None for (r, s) in zip(regions, subsets)]
    #null = [coordinates[120, 200]=labels['1000'] if r is not None else None for (r, s) in zip(regions, complete_subets)]
    
    # Position labels
    if set_labels is not None:
        # There are two situations, when set C is not on the same line with sets A and B, and when the three are on the same line.
        if abs(centers[2][1] - centers[0][1]) > tol:
            # Three circles NOT on the same line
            label_positions = [centers[0] + np.array([-radii[0] / 2, radii[0]]),
                               centers[1] + np.array([radii[1] / 2, radii[1]]),
                               centers[2] + np.array([0.0, -radii[2] * 1.1])]
            labels = [ax.text(pos[0], pos[1], txt, size=20) for (pos, txt) in zip(label_positions, set_labels)]
            labels[0].set_horizontalalignment('right')
            labels[1].set_horizontalalignment('left')
            labels[2].set_verticalalignment('top')
            labels[2].set_horizontalalignment('center')
        else:
            padding = np.mean([r * 0.1 for r in radii])
            # Three circles on the same line
            label_positions = [centers[0] + np.array([0.0, - radii[0] - padding]),
                               centers[1] + np.array([0.0, - radii[1] - padding]),
                               centers[2] + np.array([0.0, - radii[2] - padding])]
            labels = [ax.text(pos[0], pos[1], txt, size='large', ha='center', va='top') for (pos, txt) in zip(label_positions, set_labels)]
    else:
        labels = None
    return Venn3(patches, subset_labels, labels)

#### End __venn3.py
if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    import numpy as np
    #from matplotlib_venn import venn3, venn3_circles
    plt.figure(figsize=(4,4))
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    v.get_patch_by_id('100').set_alpha(1.0)
    v.get_patch_by_id('100').set_color('white')
    v.get_label_by_id('100').set_text('Unknown')
    v.get_label_by_id('A').set_text('Set "A"')
    c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    c[0].set_lw(1.0)
    c[0].set_ls('dotted')
    plt.title("Sample Venn diagram")
    plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
                ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    
    plt.show()
