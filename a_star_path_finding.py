import heapq

class Cell(object):
    def __init__(self, x, y):
        """
        Initialize new cell

        @param x cell x coordinate
        @param y cell y coordinate
        """
        self.x = x
        self.y = y
        self.parent = None
        self.g = 0
        self.h = 0
        self.f = 0

class AStar(object):
    def __init__(self):
        self.opened = []
        heapq.heapify(self.opened)
        self.closed = set()
        self.cells = []
        self.grid_height = 6
        self.grid_width = 6

    def init_grid(self):
        basic = ((0, 5), (1, 0), (1, 1), (1, 5), (2, 3), 
                 (3, 1), (3, 2), (3, 5), (4, 1), (4, 4), (5, 1))
        for x in range(self.grid_width):
            for y in range(self.grid_height):
                if (x, y) in basic:
                    self.cells.append(Cell(x, y))
        self.start=self.cells[0]
        self.end=self.cells[-1]

    def get_heuristic(self, cell1, cell2):
        """
        Compute the heuristic value H for a cell: distance between
        this cell and the ending cell multiply by 10.

        @param cell
        @returns heuristic value H
        """
        return numpy.sqrt((cell1.x - cell2.x)**2 + (cell1.y - cell2.y)**2)

    def get_cell(self, x, y):
        """
        Returns a cell from the cells list

        @param x cell x coordinate
        @param y cell y coordinate
        @returns cell
        """
        return self.cells[x * self.grid_height + y]

    def get_adjacent_cells(self, cell):
        """
        Returns adjacent cells to a cell. Clockwise starting
        from the one on the right.

        @param cell get adjacent cells for this cell
        @returns adjacent cells list 
        """
        cells = []
        if cell.x < self.grid_width-1:
            cells.append(self.get_cell(cell.x+1, cell.y))
        if cell.y > 0:
            cells.append(self.get_cell(cell.x, cell.y-1))
        if cell.x > 0:
            cells.append(self.get_cell(cell.x-1, cell.y))
        if cell.y < self.grid_height-1:
            cells.append(self.get_cell(cell.x, cell.y+1))
        return cells

    def display_path(self):
        cell = self.end
        while cell.parent is not self.start:
            cell = cell.parent
            print 'path: cell: %d,%d' % (cell.x, cell.y)

    def compare(self, cell1, cell2):
        """
        Compare 2 cells F values

        @param cell1 1st cell
        @param cell2 2nd cell
        @returns -1, 0 or 1 if lower, equal or greater
        """
        if cell1.f < cell2.f:
            return -1
        elif cell1.f > cell2.f:
            return 1
        return 0
    
    def update_cell(self, adj, cell):
        """
        Update adjacent cell

        @param adj adjacent cell to current cell
        @param cell current cell being processed
        """
        adj.g = cell.g*-1
        adj.h = self.get_heuristic(adj, cell)
        adj.parent = cell
        adj.f = -1*(adj.h + adj.g)

    def process(self):
        # add starting cell to open heap queue
        heapq.heappush(self.opened, (self.start.f, self.start))
        while len(self.opened):
            # pop cell from heap queue 
            f, cell = heapq.heappop(self.opened)
            # add cell to closed list so we don't process it twice
            self.closed.add(cell)
            # if ending cell, display found path
            #check distance to all other cells
            # get adjacent cells for cell
            for (adj_cell.f, adj_cell) in self.opened:
                self.update_cell(adj_cell, cell)
                if adj_cell.g < cell.g:
                    # add adj cell to open list
                    heapq.heappush(self.opened, (adj_cell.f, adj_cell))
                else:
                    pass
        self.display_path()

a = AStar()
a.init_grid()
a.process()

