# Shubert-Piyavski-Method for approximating the location of the global maximum(or minimum) inside a closed interval

# function to be analyzed
f = lambda z: -z ** 4 + 4 * z ** 3 + 30 * z ** 2 - 50 * z + 200
# Libshitz constant
slope = 450
eps = 0.001


# points inside the initial boundaries [a,b], such that x1 < x2
# the function returns a tuple of the cross point between the lines (x,y):
# f(x1)+slope*(x-x1) and f(x2)-slope*(x-x2), also returns x1 and x2
def find_collision_point(x1, x2):
    if x1 > x2:
        x1, x2 = x2, x1
    y = (f(x1) + f(x2) + slope * (x2 - x1)) / 2
    x = (y - f(x1)) / slope + x1
    return x, y, x1, x2


# finds if the error is insignificant
def is_not_good_approximation(expected, real):
    return abs(expected - real) > eps


# finds an approximation of the global minimum's x value
def shubert_piyavski_method(a, b):
    start = (a + b) / 2
    vertices = [find_collision_point(a, start), find_collision_point(start, b)]
    vertices.sort(key=lambda vertex: f(vertex[0]), reverse=True)

    while is_not_good_approximation(f(vertices[0][0]), vertices[0][1]):
        max_y_vertex = vertices.pop(0)
        vertices.append(find_collision_point(max_y_vertex[2], max_y_vertex[0]))
        vertices.append(find_collision_point(max_y_vertex[0], max_y_vertex[3]))
        vertices.sort(key=lambda tup: f(tup[0]), reverse=True)
    return vertices.pop()


print(shubert_piyavski_method(-5, 7))
