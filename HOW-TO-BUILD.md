Building the te
===============

```bash
cd <project_root>
mkdir -p build
cd build
cmake ..
cmake --build . --target all -- -j2 (amount of cpu's to use)
```
------------------- OR (e.g.)
```bash
mkdir -p build/CodeBlocks
cd build/CodeBlocks
cmake ../../ -G "CodeBlocks - Unix Makefiles"
codeblocks odb.cbp
```
------------------- OR

use a CMAKE Generator suitable to you environment


./synthese > ~/OpenSCAD/Getriebe/3LagenSynthese.scad
