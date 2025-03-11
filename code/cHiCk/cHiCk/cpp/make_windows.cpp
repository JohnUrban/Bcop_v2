#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

void makewindows(const std::string& BED, const std::string& G, int W, int S) {
    std::ifstream bedFile(BED);
    std::ifstream gFile(G);
    
    if (!bedFile.is_open() || !gFile.is_open()) {
        std::cerr << "Error opening files." << std::endl;
        return;
    }

    // Get Bed coords
    std::string bedLine;
    std::getline(bedFile, bedLine);
    std::istringstream bedStream(bedLine);
    
    std::vector<std::string> bedCoords;
    std::string coord;
    while (bedStream >> coord) {
        bedCoords.push_back(coord);
    }

    std::string CHR = bedCoords[0];
    int START = std::stoi(bedCoords[1]);
    int END = std::stoi(bedCoords[2]);

    // Get G
    std::string gLine;
    while (std::getline(gFile, gLine)) {
        std::istringstream gStream(gLine);
        std::vector<std::string> gCoords;
        while (gStream >> coord) {
            gCoords.push_back(coord);
        }

        if (gCoords[0] == CHR) {
            int GEND = std::stoi(gCoords[1]);

            // Get coords backwards from START
            int NEWSTART = START - S ;
            int NEWEND = NEWSTART + W;
            while (NEWSTART >= 0) {
                std::cout << CHR << "\t" << NEWSTART << "\t" << NEWEND << std::endl;
                NEWSTART -= S;
                NEWEND = NEWSTART + W;
            }

            // Get coords forwards from END
            NEWSTART = END;
            NEWEND = NEWSTART + W;
            while (NEWEND <= GEND) {
                std::cout << CHR << "\t" << NEWSTART << "\t" << NEWEND << std::endl;
                NEWSTART += S;
                NEWEND = NEWSTART + W;
            }

            break; // Stop searching for the chromosome once found
        }
    }

    bedFile.close();
    gFile.close();
}


int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <BED file> <G file> <Window size> <Step size>" << std::endl;
        return 1;
    }

    std::string bedFile = argv[1];
    std::string gFile = argv[2];
    int windowSize = std::stoi(argv[3]);
    int stepSize = std::stoi(argv[4]);

    makewindows(bedFile, gFile, windowSize, stepSize);

    return 0;
}
