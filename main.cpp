#pragma once

#include "src/lab1/lab1.hpp"
#include "src/lab2/lab2.hpp"
#include "src/lab3/lab3.hpp"
/* in dev
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
    flag = true;

    while (true) {
        std::cout << "0 - Exit" << std::endl;
        std::cout << "1. - Task 1: " << std::endl;
        std::cout << "2. - Task 2: " << std::endl;
        std::cout << "3. - Task 3: " << std::endl;
        std::cout << "4. - Task 4: " << std::endl;
        std::cout << "5. - Task 5: " << std::endl;
        std::cout << "6. - Task 6: " << std::endl;
        std::cout << "7. - Task 7: " << std::endl;
        std::cout << "8. - Task 8: " << std::endl;
        std::cout << "9. - Task 9: " << std::endl;
        std::cout << "10. - Task 10: " << std::endl;
        std::cout << "11. - Task 11: " << std::endl;
        std::cout << "12. - Task 12: " << std::endl;
        std::cout << "13. - Task 13: " << std::endl;
        std::cout << "14. - Task 14: " << std::endl;
        std::cout << "15. - Task 15: " << std::endl;
        std::cout << "16. - Task 16: " << std::endl;
        std::cout << "17. - Task 17: " << std::endl;
        std::cout << "18. - Task 18: " << std::endl;
        std::cout << ">>";

        std::cin >> switcher;

        switch(switcher) {
            case 0: {
                flag = !flag;
            }
            case 1: {
                lab1_start();
                break;
            }
            case 2: {
                lab2_start();
                break;
            }
            case 3: {
                lab3_start();
                break;
            }/* in dev
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
                std::cout << "No such option, try again!" << std::endl;
                break;
            }
        }
    }
    return 0;
}