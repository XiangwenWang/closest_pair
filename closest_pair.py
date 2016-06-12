"""
A program designed to compute the closest pair in d-dimensional space.
The algorithm efficiency is between O(n) to O(log n).
Author: Xiangwen Wang
"""
import random


class closest_pair:
    _break_condition = 0.02

    def __init__(self, pts):
        if type(pts) is set:
            self.pts = pts
        else:
            self.pts = set(pts)
        self.pts_copy = self.pts.copy()
        self.n = len(self.pts)
        self.listDuplicated = False
        if self.n < len(pts):
            return self.findDuplicate(pts)
        self.dim = len(random.choice(tuple(self.pts)))
        self.origin = self.getOrigin()
        # print self.dim, self.origin

    def findDuplicate(self, pts):
        """
        If there are duplicate points, then they are the closest pair
        """
        for x in pts:
            if x in self.pts:
                self.pts.remove(x)
            else:
                break
        self.listDuplicated = (True, [0.0, [x, x]])

    def distance(self, x, y):
        return sum(map(lambda p, q: (p - q) ** 2, x, y)) ** 0.5

    def getOrigin(self):
        dim = self.dim
        origin = [False] * self.dim
        for x in self.pts:
            for i in xrange(dim):
                if origin[i] is False or x[i] < origin[i]:
                    origin[i] = x[i]
        return origin

    def calGridSize(self):
        """
        Determine the Grid Size b, which is the cloest distance
        from a random chosen point to the rest points
        """
        x = random.choice(tuple(self.pts))
        close_di = min([self.distance(x, xi) for xi in self.pts if xi != x])
        # print close_di
        return close_di

    def removePts(self):
        """
        If one point has an empty neighboor, remove the point
        It does not belong to the cloeset pair
        """
        self.gridsize = self.calGridSize()
        # print self.gridsize
        self.mesh()
        for x in self.pts.copy():
            coord_x = self.pts2grid[x]
            if len(self.grid2pts[coord_x]) == 1:
                if self.emptyNeighb(coord_x)[0]:
                    self.pts.remove(x)
                    # print self.emptyNeighb(x)[0],len(self.pts)
        self.n = len(self.pts)

    def emptyNeighb(self, coord):
        """
        Determine whether the given point has an empty neighboor
        """
        neighboors = set()
        neighboors.add(coord)
        for i in xrange(self.dim):
            for neigh_coord in neighboors.copy():
                for j in (-1, 1):
                    new_neigh_coord = list(neigh_coord)
                    new_neigh_coord[i] += j
                    new_neigh_coord = tuple(new_neigh_coord)
                    if new_neigh_coord not in neighboors:
                        neighboors.add(new_neigh_coord)
        neighboors.remove(coord)
        for y in neighboors.copy():
            if -1 in y:
                neighboors.remove(y)
        for y in neighboors:
            if y in self.grid2pts:
                return (False, neighboors)
        return (True, neighboors)

    def mesh(self):
        """
        Meshing the points using the gridsize calculated above
        """
        pts2grid = {}
        grid2pts = {}
        b = self.gridsize
        ori = self.origin
        for x in self.pts:
            coord = tuple(map(lambda y, z: int((y - z) / b), x, ori))
            pts2grid[x] = coord
            if coord in grid2pts:
                grid2pts[coord].append(x)
            else:
                grid2pts[coord] = [x]
        self.grid2pts = grid2pts
        self.pts2grid = pts2grid
        # print self.grid2pts
        # print self.pts2grid

    def residClosest(self):
        """
        After several meshes and removing the empty-neighboor points
        find the cloest pair in the remained points
        """
        if self.n == 2:
            pt_pair = list(self.pts)
            dist = self.distance(pt_pair[0], pt_pair[1])
            return [dist, pt_pair]
        # return self.BFclosest_d()
        closest_d = False
        closestPair = False
        for x in self.pts:
            neigh_pts = set()
            coord_x = self.pts2grid[x]
            for xi in self.grid2pts[coord_x]:
                if xi != x:
                    neigh_pts.add(xi)
            for coord in self.emptyNeighb(coord_x)[1]:
                if coord in self.grid2pts:
                    for xi in self.grid2pts[coord]:
                        neigh_pts.add(xi)
            for xi in neigh_pts:
                dist = self.distance(xi, x)
                if closest_d is False or dist < closest_d:
                    closestPair = [x, xi]
                    closest_d = dist
        return [closest_d, closestPair]

    def BFclosest_d(self):
        """
        Finding the cloest pair using the brute force method
        """
        if self.listDuplicated:
            return self.listDuplicated[1]
        pts_exists = set()
        closest_d = False
        closestPair = False
        for i in self.pts:
            for j in pts_exists:
                dist = self.distance(i, j)
                if closest_d is False or dist < closest_d:
                    closestPair = [i, j]
                    closest_d = dist
            pts_exists.add(i)
        self.bf_d = [closest_d, closestPair]
        return [closest_d, closestPair]

    def neighbSear(self):
        """
        Finding the cloest pair using the mesh-neighboor-search method
        """
        if self.listDuplicated:
            return self.listDuplicated[1]
        n = self.n
        while n > 2:
            self.removePts()
            if n - self.n <= self._break_condition * n:
                break
            # print grid.n, n
            n = self.n
        self.d = self.residClosest()
        return self.d

    @staticmethod
    def genRndPts(n, ori, scale):
        """
        Generate random points for testing
        """
        import sys
        if len(sys.argv) == 2:
            seed = int(sys.argv[1])
            random.seed(seed)
        PtsSet = set()
        for i in xrange(n):
            pt = tuple(map(lambda p, q: p + random.random() * q, ori, scale))
            PtsSet.add(pt)
        return PtsSet

    @staticmethod
    def test():
        import time
        pts = closest_pair.genRndPts(2000, [101.6, 26.75, -12.56, 1992.5],
                                     [200, 300, 100, 413])
        grid = closest_pair(pts.copy())
        begin = time.time()
        d = grid.BFclosest_d()[0]
        end = time.time()
        print -1, d, end - begin, 'bf-solution'
        for i in xrange(100):
            grid = closest_pair(pts.copy())
            begin = time.time()
            d = grid.neighbSear()[0]
            end = time.time()
            print i, d, end - begin

    @staticmethod
    def test_list():
        import time
        pts = list(closest_pair.genRndPts(2000, [101.6, 26.75, -12.56, 1992.5],
                                          [200, 300, 100, 413]))
        # pts = pts + pts[0:1]
        grid = closest_pair(pts[:])
        begin = time.time()
        d = grid.BFclosest_d()[0]
        end = time.time()
        print -1, d, end - begin, 'bf-solution'
        for i in xrange(1000):
            grid = closest_pair(pts[:])
            begin = time.time()
            d = grid.neighbSear()[0]
            end = time.time()
            print i, d, end - begin


def main():
    # closest_pair.test()
    closest_pair.test_list()


if __name__ == '__main__':
    main()
