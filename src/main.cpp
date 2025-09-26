#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include "OccupancyGrid.h"
#include "PathPlannerInterface.h"


/**
 * @class AStarPlanner
 * @brief Implements the PathPlannerInterface using the A* algorithm.
 */
class AStarPlanner : public PathPlannerInterface {
public:
    Trajectory getCollisionFreePath(const OccupancyGrid& grid,
                                      const Eigen::Vector3f& start,
                                      const Eigen::Vector3f& end) override;

private:
    /**
    * @struct Node
    * @brief Represents a node in the A* search space.
    */
    struct Node {
        int y;
        int x;
        double h;
        double g;
        double f;
        Node* parent;

        Node(int y, int x, double h, double g, Node* parent = nullptr) :
            
        y(y), x(x), h(h), g(g), f(g + h), parent(parent) {}
    };

    /**
    * @struct Node
    * @brief To make the priority_queue open_list a min-heap
    */
    struct CompareNodeCosts {
        bool operator()(const Node* a, const Node* b) const {
            return a->f > b->f;
        }
    };

    /**
     * @brief Creates a new grid with obstacles inflated by the robot's radius.
     * @param originalGrid The initial map.
     * @param robotRadius The robot's radius.
     * @return A new OccupancyGrid with expanded obstacles.
     */
    static OccupancyGrid createInflatedGrid(const OccupancyGrid& originalGrid, float robotRadius);

    /**
     * @brief Checks if a cell is valid (within grid bounds and not an obstacle).
     * @param y The row index of the cell.
     * @param x The column index of the cell.
     * @param grid The occupancy grid (should be the inflated grid).
     * @return True if the cell is valid, false otherwise.
     */
    static bool is_valid(int y, int x, const OccupancyGrid& grid);

    /**
     * @brief Calculates the Euclidean distance heuristic.
     * @return The heuristic distance.
     */
    static double heuristic(int y1, int x1, int y2, int x2);
};

// A* Implementation
Trajectory AStarPlanner::getCollisionFreePath(const OccupancyGrid& grid,
                                              const Eigen::Vector3f& start,
                                              const Eigen::Vector3f& end) {

    const float robotRadius = 0.3;
    OccupancyGrid inflatedGrid = createInflatedGrid(grid, robotRadius);

    int start_y = static_cast<int>(start.y() / inflatedGrid.resolution_m);
    int start_x = static_cast<int>(start.x() / inflatedGrid.resolution_m);
    int end_y = static_cast<int>(end.y() / inflatedGrid.resolution_m);
    int end_x = static_cast<int>(end.x() / inflatedGrid.resolution_m);

    Node* start_node = new Node(start_y, start_x, heuristic(start_y, start_x, end_y, end_x), 0);
    
    std::priority_queue<Node*, std::vector<Node*>, CompareNodeCosts> open_list;    
    std::unordered_map<int, Node*> node_map;
    open_list.push(start_node);
    node_map[inflatedGrid.get1DIndex(start_x, start_y)] = start_node;

    Trajectory path;

    while (!open_list.empty()) {
        Node* current_node = open_list.top();
        open_list.pop();

        if (current_node->y == end_y && current_node->x == end_x) {
            Node* temp = current_node;
            while (temp != nullptr) {
                path.push_back(Eigen::Vector3f(temp->x * inflatedGrid.resolution_m, temp->y * inflatedGrid.resolution_m, 0));
                temp = temp->parent;
            }
            std::reverse(path.begin(), path.end());
            break;
        }

        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dy == 0 && dx == 0) continue;

                int new_y = current_node->y + dy;
                int new_x = current_node->x + dx;

                if (is_valid(new_y, new_x, inflatedGrid)) {
                    double new_g = current_node->g + std::sqrt(std::pow(dy, 2) + std::pow(dx, 2));
                    int key = inflatedGrid.get1DIndex(new_x, new_y);

                    if (node_map.find(key) == node_map.end() || new_g < node_map[key]->g) {
                        Node* neighbor = new Node(new_y, new_x, heuristic(new_y, new_x, end_y, end_x), new_g, current_node);
                        open_list.push(neighbor);
                        node_map[key] = neighbor;
                    }
                }
            }
        }
    }

    for (auto const& [key, val] : node_map) {
        delete val;
    }
    
    return path;
}

OccupancyGrid AStarPlanner::createInflatedGrid(const OccupancyGrid& originalGrid, float robotRadius_m) {
    OccupancyGrid inflatedGrid = originalGrid;
    int radiusInCells = static_cast<int>(std::ceil(robotRadius_m / originalGrid.resolution_m));

    for (int y = 0; y < originalGrid.height; ++y) {
        for (int x = 0; x < originalGrid.width; ++x) {
            if (originalGrid.data[originalGrid.get1DIndex(x, y)]) {
                for (int dy = -radiusInCells; dy <= radiusInCells; ++dy) {
                    for (int dx = -radiusInCells; dx <= radiusInCells; ++dx) {
                        if (std::sqrt(dx*dx + dy*dy) <= radiusInCells) {
                            int newY = y + dy;
                            int newX = x + dx;
                            if (newY >= 0 && newY < inflatedGrid.height && newX >= 0 && newX < inflatedGrid.width) {
                                inflatedGrid.data[newY * inflatedGrid.width + newX] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    return inflatedGrid;
}

bool AStarPlanner::is_valid(int y, int x, const OccupancyGrid& grid) {
    return y >= 0 && y < grid.height && x >= 0 && x < grid.width && !grid.data[grid.get1DIndex(x, y)];
}

double AStarPlanner::heuristic(int y1, int x1, int y2, int x2) {
    return std::sqrt(std::pow(y1 - y2, 2) + std::pow(x1 - x2, 2));
}

int main() {
    OccupancyGrid grid{};
    // Sample wall obstacle
    for(int i = 50; i < 250; ++i) {
        grid.setAsObstacle(i * grid.resolution_m, 150 * grid.resolution_m);
    }

    std::unique_ptr<PathPlannerInterface> planner = std::make_unique<AStarPlanner>();

    auto start = Eigen::Vector3f(2.0, 2.0, 0.0);
    auto end = Eigen::Vector3f(14.0, 14.0, 0.0);

    Trajectory path = planner->getCollisionFreePath(grid, start, end);

    if (path.empty()) {
        std::cout << "No path found!" << std::endl;
    } else {
        std::cout << "Path found! Path:" << std::endl;
        for (const auto& point : path) {
            std::cout << "  (" << point.x() << ", " << point.y() << ")" << std::endl;
        }
    }

    return 0;
}