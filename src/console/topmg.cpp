//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "console/topmg_argument.hpp"
#include "console/topmg_process.hpp"

int main(int argc, char* argv[]) {
  // toppic::log_level = 2;
  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::cout << "The starting time of the program: " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  std::cout << std::setprecision(10);
  toppic::Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);
  if (!success) {
    return 1;
  }
  std::map<std::string, std::string> arguments = argu_processor.getArguments();
  
  std::vector<std::string> spec_file_lst = argu_processor.getSpecFileList();

  std::sort(spec_file_lst.begin(), spec_file_lst.end());

  return toppic::TopMGProgress_multi_file(arguments, spec_file_lst);
}
