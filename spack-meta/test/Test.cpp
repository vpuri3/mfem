#include <gtest/gtest.h>
#include "AStyleTests.cpp"

int main(int argc, char **argv)
{
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}