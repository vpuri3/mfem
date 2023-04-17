#include <iostream>

using namespace std;

int main()
{
    int a  = 5;
    int *p = &a;

    cout << a  << endl; // 5
    cout << p  << endl; // 0x7fff0d4c2934
    cout << &a << endl; // 0x7fff0d4c2934
    cout << *p << endl; // 5

}
