/*Klems Annual Series on Output, employment, and capital*/

use "C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\replication_elasticities\klem.dta" , clear

cd "C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\code_to_eric\"
keep year qout qktemp qltemp qm ind po

sort ind year
gen ind2=1 if ind==1
replace ind2=2 if ind==2 | ind==3 | ind==5
replace ind2=3 if ind==4
replace ind2=4 if ind==6
replace ind2=5 if ind==7 | ind==8
replace ind2=6 if ind==9
replace ind2=7 if ind==10 | ind==18
replace ind2=ind-3 if ind>=11 & ind<=17
replace ind2=ind-4 if ind>=19 & ind<=29
replace ind2=26 if ind==30 | ind==31
replace ind2=ind-5 if ind>=32 & ind<=35

collapse (sum) qout qktemp qltemp qm (mean) po, by(ind2 year)
duplicates drop
ren ind2 ind

preserve
keep if year==1997
g sales=po*qout
egen tsales=sum(sales)
g shry=sales/tsales
egen ttsales=sum(shry)
keep shry
outsheet using output_shares_klems.csv, replace delim(,) nonames
restore

sort year ind
reshape wide qout qktemp qltemp qm, i(year) j(ind)

preserve	  
sort year
keep qout*
outsheet using data_output_klems.csv, replace delim(,) nonames
restore

preserve	  
sort year
keep qktemp*
outsheet using data_capital_klems.csv, replace delim(,) nonames
restore

preserve	  
sort year
keep qltemp*
outsheet using data_labor_klems.csv, replace delim(,) nonames
restore

preserve	  
sort year
keep qm*
outsheet using data_material_klems.csv, replace delim(,) nonames
restore

