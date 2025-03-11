#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

void grablinks(const std::string& CHR, int A, int B, const std::string& BEDPE) {
    std::ifstream bedpeFile(BEDPE);
    
    if (!bedpeFile.is_open()) {
        std::cerr << "Error opening file: " << BEDPE << std::endl;
        return;
    }

    std::string line;
    while (std::getline(bedpeFile, line)) {
        std::istringstream bedpeStream(line);
        std::vector<std::string> bedpeFields;
        std::string field;
        while (bedpeStream >> field) {
            bedpeFields.push_back(field);
        }

        // Check conditions and print if satisfied
        if (bedpeFields[0] == CHR && bedpeFields[3] == CHR &&
            ((std::stoi(bedpeFields[1]) >= A && std::stoi(bedpeFields[1]) <= B) ||
             (std::stoi(bedpeFields[2]) >= A && std::stoi(bedpeFields[2]) <= B) ||
             (std::stoi(bedpeFields[4]) >= A && std::stoi(bedpeFields[4]) <= B) ||
             (std::stoi(bedpeFields[5]) >= A && std::stoi(bedpeFields[5]) <= B))) {
            std::cout << line << std::endl;
        }
    }

    bedpeFile.close();
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <CHR> <A> <B> <BEDPE file>" << std::endl;
        return 1;
    }

    std::string CHR = argv[1];
    int A = std::stoi(argv[2]);
    int B = std::stoi(argv[3]);
    std::string BEDPE = argv[4];

    grablinks(CHR, A, B, BEDPE);

    return 0;
}
