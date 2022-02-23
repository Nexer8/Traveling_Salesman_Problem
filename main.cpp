#include "Algorithm.hpp"
#include "utils/GUI.hpp"

using namespace std;

int main() {
    cout << "****************************************" << endl;
    cout << "Mariusz Wiśniewski" << endl;

    int rc = 0;
    while (rc != EXIT) {
        rc = GUI::run();
    }

    return 0;
}