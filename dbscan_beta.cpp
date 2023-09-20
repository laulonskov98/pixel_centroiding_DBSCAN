#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath> 
#include <chrono> //for measuring time
#include <iomanip> //for setprecision
#include <cstdio>


//struct for a point
struct Point {
    double shot; //shot id
    double x;
    double y;
    double z; //tof transformed to z coordinate
    double tof;
    double tot;
    int label; // Cluster label
    bool visited; // Flag to mark if the point has been visited during clustering
};

//struct for a cluster
struct Cluster {
    int id;          // Cluster ID
    double shot;     // Shot ID
    double avgX;     // Average x coordinate
    double avgY;     // Average y coordinate
    double avgTOF;   // Average tof
    double maxTOT;   // Maximum tot
    std::vector<Point> points;  // Points in the cluster
};



double interpolate(const std::vector<double>& xVals, const std::vector<double>& yVals, double x) {
    int i = 0;

    if (x >= xVals[xVals.size() - 2]) {
        i = xVals.size() - 2; //edge case for last element
    } else {
        // Find index i so that xVals[i] <= x < xVals[i+1]
        while (x > xVals[i + 1]) {
            i++;
        }
    }   

    //initialize variables for linear interpolation
    double xLow = xVals[i], yLow = yVals[i];
    double xHigh = xVals[i + 1], yHigh = yVals[i + 1];

    //throw error if x is outside of the range since no extrapolation handled currently
    if (x < xLow || x > xHigh) {
        std::cerr << "Error: x value outside of range for interpolation: " << x << std::endl;
        return 0;
    }

    double dydx = (yHigh - yLow) / (xHigh - xLow); // compute slope

    return yLow + dydx * (x - xLow); // Linear interpolation
}


// Calculate Euclidean distance between two points
double calculateDistance(const Point& p1, const Point& p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                (p1.y - p2.y) * (p1.y - p2.y) +
                (p1.z - p2.z) * (p1.z - p2.z));
}

// Find neighboring points within the epsilon distance
std::vector<int> findNeighbors(const std::vector<Point>& points, int pointIndex, double epsilon) {
    std::vector<int> neighbors;
    for (int i = 0; i < points.size(); ++i) {
        if (i == pointIndex) continue;
        double distance = calculateDistance(points[pointIndex], points[i]);
        if (distance <= epsilon) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

// Expand a cluster starting from a seed point
void expandCluster(std::vector<Point>& points, int pointIndex, int cluster, double epsilon, int minPoints) {
    std::vector<int> seeds = findNeighbors(points, pointIndex, epsilon);
    if (seeds.size() < minPoints) {
        points[pointIndex].label = -1; // Mark as noise. Should already be marked by default, but just to be safe.
        return;
    }
    points[pointIndex].label = cluster;
    points[pointIndex].visited = true;

    
    for (int i = 0; i < seeds.size(); ++i) {
        int seedIndex = seeds[i];
        if (!points[seedIndex].visited) {
            expandCluster(points, seedIndex, cluster, epsilon, minPoints);
        }
    }
}


//dbscan algorithm
void dbscan(std::vector<Point>& points, double epsilon, int minPoints, int& cluster) {
    for (int i = 0; i < points.size(); ++i) {
        if (!points[i].visited) {
            expandCluster(points, i, cluster, epsilon, minPoints);
            if (points[i].label != -1) {
                cluster++;
            }
        }
    }
}



//########################################################################################
//########################################################################################
//########################################################################################

// Main function
int main(int argc, char* argv[]) {
    if (argc != 3 & argc != 4) {
        std::cerr << "Usage: " << argv[0] << "<input_csv_file> <output_csv_file> <optional correctionfile>" << std::endl;
        return 1;
    }

    const char* inputFileName = argv[1]; // Get the input file name from the command line
    const char* outputFileName = argv[2]; // Get the output file name from the command line
    const char* correctionFile = NULL; 
    if (argc == 4) {
        correctionFile = argv[3]; // Get the correction file name from the command line
    }


    //DBSCAN parameters
    const double EPSILON = 2.0; // Epsilon neighborhood distance threshold
    const int MINPOINTS = 1;    // Min number of points required to form a cluster (SELF EXCLUDED)
    
    //tof parameters
    const double TOFTHRESHOLD = 2e-4; //tof threshold for filtering points
    const double TOFTTRANSFORM = 81920*(25./4096)*1E-9; //weird chemistry science stuff magic number maybe??
    const double TOFTTOZ = EPSILON / TOFTTRANSFORM; //tof to z coordinate transformation factor

    const int DOUBLEPRECISSION = 17; //number of digits for double precission
    const int EXPECTEDNUMBERPOINTS = 5e5; //expected number of datapoints in input file

    
    //linear interpolation of correction file (if provided)
    std::vector<double> tofVals;
    std::vector<double> correctionVals;
    if (correctionFile != NULL) {
        std::ifstream file(correctionFile);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open the file: " << correctionFile << std::endl;
            return 1;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double tof;
            double correction;
            char delim; // Variable to store the delimiter (,) character
            //read correction file assuming format: tof,correction and no header
            if (iss >> tof >> delim >> correction) {
                tofVals.push_back(tof);
                correctionVals.push_back(correction);
            } else {
                std::cerr << "Error: Invalid data format in line: " << line << std::endl;
            }
        }
    }


    // Open the input file C style
    FILE* file = fopen(inputFileName, "r");
    if (!file) {
        std::cerr << "Error: Could not open the file: " << inputFileName << std::endl;
        return 1;
    }
    
    // Create a vector to accumulate the points and reserve space to avoid reallocation
    std::vector<Point> points;
    points.reserve(EXPECTEDNUMBERPOINTS);

    // Read the header line C style
    char header[256]; 
    if (fgets(header, sizeof(header), file) == NULL) {
        std::cerr << "Error: Failed to read the header" << std::endl;
        fclose(file);
        return 1;
    }

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); //measure time for reading input file

    int shotOffset; //Variable to store the shot offset from the input file
    bool firstShot = true; //Variable to store the shot offset from the input file

    // Read data using fscanf
    while (!feof(file)) {
        Point point;
        //read input file assuming format: shot,x,y,tof,tot and no header (C style)
        if (fscanf(file, "%*d,%lf,%lf,%lf,%lf,%lf", &point.shot, &point.x, &point.y, &point.tof, &point.tot) == 5) {
            //Create new shot ids starting from 1 (instead of the numbers from the input file)
            if (point.tof > TOFTHRESHOLD) {continue;} //skip points with tof above threshold
            else {point.z = point.tof*TOFTTOZ;} //transform tof to z coordinate

            //Create new shot ids starting from 1 (instead of the numbers from the input file)
            if (firstShot) {
                shotOffset = point.shot;
                firstShot = false;
            }
            point.shot = point.shot-shotOffset+1; //shift shot ids to start from 1

            // Initialize point label to -1 (unassigned)
            point.label = -1;
            point.visited = false;


            //linear interpolation of correction file (if provided)
            if (correctionFile != NULL) {
                double tofCorrection = interpolate(tofVals, correctionVals, point.tot);
                point.tof = point.tof - tofCorrection;
            }
            points.push_back(point);

        } else {
            std::cerr << "Error: Invalid data format in line." << std::endl;
        }
    }



    auto finish = std::chrono::high_resolution_clock::now(); //measure time for reading input file
    std::chrono::duration<double> elapsed = finish - start; //measure time for reading input file
    std::cout << "Elapsed time for reading and processing input file: " << elapsed.count()*1000 << " ms\n"; //measure time for reading input file


    start = std::chrono::high_resolution_clock::now(); //measure time for grouping data
    std::map<double, std::vector<Point>> shotToPoints;

    // Group points by 'shot'
    for (const Point& point : points) {
        shotToPoints[point.shot].push_back(point);
    }

    finish = std::chrono::high_resolution_clock::now(); //measure time for grouping data
    elapsed = finish - start; //measure time for grouping data
    std::cout << "Elapsed time for grouping points: " << elapsed.count()*1000 << " ms\n"; //measure time for grouping data

    start = std::chrono::high_resolution_clock::now(); //measure time for DBSCAN algorithm

    // Perform DBSCAN on each shot
    int clusterId = 0;
    for (auto& pair : shotToPoints) {
        dbscan(pair.second, EPSILON, MINPOINTS,clusterId);
    }

    finish = std::chrono::high_resolution_clock::now(); //measure time for DBSCAN algorithm
    elapsed = finish - start; //measure time for DBSCAN algorithm
    std::cout << "Elapsed time for DBSCAN: " << elapsed.count()*1000 << " ms\n"; //measure time for DBSCAN algorithm



    start = std::chrono::high_resolution_clock::now(); //measure time for cluster processing

    // Create a map to accumulate the clusters
    std::map<int, Cluster> clusters; // Map cluster ID to Cluster structure
    // Loop through the points and assign them to clusters based on their labels
    for (const auto& pair : shotToPoints) {
        double shot = pair.first;
        for (const Point& point : pair.second) {
            if (point.label != -1) { // Ignore noise points
                int clusterLabel = point.label;
                // Check if the cluster exists in the map, and if not, create it
                if (clusters.find(clusterLabel) == clusters.end()) {
                    clusters[clusterLabel].id = clusterLabel; // Initialize cluster ID
                    clusters[clusterLabel].shot = shot; // Initialize shot
                    clusters[clusterLabel].maxTOT = point.tot; // Initialize maxTOT
                }
                clusters[point.label].points.push_back(point); // Add point to cluster
                // Update the maximum tot value for the cluster
                if (point.tot > clusters[point.label].maxTOT) {
                    clusters[point.label].maxTOT = point.tot;
                }
            }
        }
    }
    
    // Calculate averages and for each cluster
    for (auto& pair : clusters) {
        Cluster& cluster = pair.second;
        double sumX = 0.0, sumY = 0.0, sumTOF = 0.0, sumWeight = 0.0; //initialize variables for weighted average

        if (cluster.points.empty()) {continue;} //skip empty clusters

        for (const Point& point : cluster.points) {
            sumX += point.x * point.tot;
            sumY += point.y * point.tot;
            sumTOF += point.tof * point.tot;
            sumWeight += point.tot; //needed for weighted average
        }
        // Calculate weighted averages
        cluster.avgX = sumX / sumWeight;
        cluster.avgY = sumY / sumWeight;
        cluster.avgTOF = sumTOF / sumWeight;

    }
    finish = std::chrono::high_resolution_clock::now(); //measure time for cluster processing
    elapsed = finish - start; //measure time for cluster processing


    std::cout << "Elapsed time for cluster processing: " << elapsed.count()*1000 << " ms\n"; //measure time for cluster processing



    start = std::chrono::high_resolution_clock::now(); //measure time for writing output file
    // Create a vector to accumulate the output
    std::vector<std::string> outputData;
    //output header line
    outputData.push_back("Index,shot,x,y,tof,tot\n");
    for (auto& pair: clusters) {
        Cluster cluster = pair.second;
        std::ostringstream oss;
        oss << pair.first << "," 
            << cluster.shot << ","
            << std::setprecision(DOUBLEPRECISSION) << cluster.avgX << "," 
            << std::setprecision(DOUBLEPRECISSION) << cluster.avgY << "," 
            << std::setprecision(DOUBLEPRECISSION)<< cluster.avgTOF << ","
            << std::setprecision(DOUBLEPRECISSION)  << cluster.maxTOT << "\n";

        outputData.push_back(oss.str());

    }

    // Open the output file for writing
    std::ofstream outputFile(outputFileName);

    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the output file: " << outputFileName << std::endl;
        return 1;
    }

    // Write all the accumulated data to the output file
    for (const std::string& line : outputData) {
        outputFile << line;
    }
    finish = std::chrono::high_resolution_clock::now(); //measure time for writing output file
    elapsed = finish - start; //measure time for writing output file
    std::cout << "Elapsed time for writing output file: " << elapsed.count()*1000 << " ms\n"; //measure time for writing output file
    

    std::cout << "Output written to " << outputFileName << std::endl;

    return 0;
}