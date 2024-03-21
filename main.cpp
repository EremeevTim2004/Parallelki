#pragma once

#include "src/lab1/lab1.hpp"
#include "src/lab2/lab2.hpp"
/* in dev
#include "src/lab3/lab3.hpp"
#include "src/lab4/lab4.hpp"
#include "src/lab5/lab5.hpp"
#include "src/lab6/lab6.hpp"
#include "src/lab7/lab7.hpp"
#include "src/lab8/lab8.hpp"
#include "src/lab9/lab9.hpp"
#include "src/lab10/lab10.hpp"
#include "src/lab11/lab11.hpp"
#include "src/lab12/lab12.hpp"
#include "src/lab13/lab13.hpp"
#include "src/lab14/lab14.hpp"
#include "src/lab15/lab15.hpp"
#include "src/lab16/lab16.hpp"
#include "src/lab17/lab17.hpp"
#include "src/lab18/lab18.hpp"
*/
using namespace std;

int main() {
    int switcher = 0;
    switch(switcher) {
        case 1: {
            lab1_start();
            break;
        }
        case 2: {
            lab2_start();
            break;
        }/* in dev
        case 3: {
            lab3_start();
            break;
        }
        case 4: {
            lab4_start();
            break;
        }
        case 5: {
            lab5_start();
            break;
        }
        case 6: {
            lab6_start();
            break;
        }
        case 7: {
            lab7_start();
            break;
        }
        case 8: {
            lab8_start();
            break;
        }
        case 9: {
            lab9_start();
            break;
        }
        case 10: {
            lab10_start();
            break;
        }
        case 11: {
            lab11_start();
            break;
        }
        case 12: {
            lab12_start();
            break;
        }
        case 13: {
            lab13_start();
            break;
        }
        case 14: {
            lab14_start();
            break;
        }
        case 15: {
            lab15_start();
            break;
        }
        case 16: {
            lab16_start();
            break;
        }
        case 17: {
            lab17_start();
            break;
        }
        case 18: {
            lab18_start();
            break;
        }*/
        default: {
            break;
        }
    }
    return 0;
}