//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_GUI_HPP
#define TSP_GUI_HPP

#include "../Algorithm.hpp"

#include<string>

using namespace std;

#define EXIT (-2)
#define REPEAT 1

class GUI {
public:
    static int run();

    static string displayMenuStart();
};

#endif //TSP_GUI_HPP