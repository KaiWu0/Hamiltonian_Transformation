#/bin/bash
#This script copies ht.f90 to PP/src and modifies PP/src/Makefile.
if grep -q "ht.x" PP/src/Makefile; then
    echo PP/src/Makefile already contains keyword \"ht.x\", do nothing.
else
    cp ht.f90 PP/src
    sed -ri '/^all\ :/ s/\\$/ht.x\ \\/' PP/src/Makefile
    sed -i '/^clean\ :.*/i ht.x : ht.o libpp.a $(MODULES)\n\t$(LD) $(LDFLAGS) -o $@ \\\n\t\tht.o libpp.a $(MODULES) $(QELIBS)\n\t- ( cd ../../bin ; ln -fs ../PP/src/$@ . )\n' PP/src/Makefile
    echo "success"
fi
