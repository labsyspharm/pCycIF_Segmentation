function printDrug(drugName)
if isnumeric(drugName)
    sprintf('Working on drug %s',num2str(drugName))
else
    sprintf('Working on drug %s',drugName)
end