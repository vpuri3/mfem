#include <gtest/gtest.h>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <filesystem>

std::string exec(const char* cmd)
{
   std::array<char, 128> buffer;
   std::string result;
   std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
   if (!pipe)
   {
      throw std::runtime_error("popen() failed!");
   }
   while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
   {
      result += buffer.data();
   }
   return result;
}

TEST(AStyleTests, AStyleSrcCPP)
{
   std::string astyle_out, cmd;
   for (auto&p: std::filesystem::recursive_directory_iterator("../src"))
   {
      if (p.is_regular_file() &&
          p.path().extension().string().find("cpp") != std::string::npos)
      {
         cmd = "astyle --dry-run --options=../config/fluids-mfem.astylerc " +
               p.path().string();
         astyle_out = exec(cmd.c_str());
         EXPECT_TRUE(astyle_out.find("Formatted") == std::string::npos) << "File: " +
                                                                        p.path().string();
      }
   }
}

TEST(AStyleTests, AStyleSrcHPP)
{
   std::string astyle_out, cmd;
   for (auto&p: std::filesystem::recursive_directory_iterator("../src"))
   {
      if (p.is_regular_file() &&
          p.path().extension().string().find("hpp") != std::string::npos)
      {
         cmd = "astyle --dry-run --options=../config/fluids-mfem.astylerc " +
               p.path().string();
         astyle_out = exec(cmd.c_str());
         EXPECT_TRUE(astyle_out.find("Formatted") == std::string::npos) << "File: " +
                                                                        p.path().string();
      }
   }
}

TEST(AStyleTests, AStyleTestCPP)
{
   std::string astyle_out, cmd;
   for (auto&p: std::filesystem::recursive_directory_iterator("../test"))
   {
      if (p.is_regular_file() &&
          p.path().extension().string().find("cpp") != std::string::npos)
      {
         cmd = "astyle --dry-run --options=../config/fluids-mfem.astylerc " +
               p.path().string();
         astyle_out = exec(cmd.c_str());
         EXPECT_TRUE(astyle_out.find("Formatted") == std::string::npos) << "File: " +
                                                                        p.path().string();
      }
   }
}

TEST(AStyleTests, AStyleTestHPP)
{
   std::string astyle_out, cmd;
   for (auto&p: std::filesystem::recursive_directory_iterator("../test"))
   {
      if (p.is_regular_file() &&
          p.path().extension().string().find("hpp") != std::string::npos)
      {
         cmd = "astyle --dry-run --options=../config/fluids-mfem.astylerc " +
               p.path().string();
         astyle_out = exec(cmd.c_str());
         EXPECT_TRUE(astyle_out.find("Formatted") == std::string::npos) << "File: " +
                                                                        p.path().string();
      }
   }
}
