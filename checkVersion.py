import re

rad_varsion_file = "src/rt/VERSION"
rad_varsion = "5.3.0"


with open(rad_varsion_file, "r") as f:
    rad_varsion = f.readline().strip()
    print(rad_varsion)


cmake_list_file = "CMakeLists.txt"
regex = "(?<=project\(radiance\sVERSION\s)\w*.\w*.\w*(?=\))"
replace_new = re.sub('[a-zA-Z\s]','', rad_varsion )
# replace_new = "5.3.1"
print(replace_new)

with open(cmake_list_file, "r") as f:
    s = f.read()
with open(cmake_list_file, 'w') as f:
    s = re.sub(regex, replace_new, s)
    f.write(s)