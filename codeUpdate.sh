#!/bin/sh
rm ./albertoCode.tar.gz
tar -czvf ./albertoCode.tar.gz ./PDAnalysis/
scp ./albertoCode* abragagn@t2-ui-12:/lustre/cmswork/abragagn/BPH/
#scp -2 -P 4444 -oNoHostAuthenticationForLocalhost=yes ./albertoCode.tar.gz  abragagn@localhost:/lustre/cmswork/abragagn/BPH/

git commit -a -m "Update"
git push
