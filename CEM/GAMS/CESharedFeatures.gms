*Michael Craig 16 Apr 2020

$include EDAndUCSharedFeatures.gms

Sets
*EXISTING GENERATORS
         nonrenewegu(egu)                existing nonrenewable generators
                 ststorageegu(storageegu)
                 ltstorageegu(storageegu)
*CANDIDATE TECHNOLOGIES FOR CONSTRUCTION
         tech                            candidate technologies for new construction
*Generators
                gentechs(tech)                                  generating techs specifically not DACS or storage
         thermaltech(tech)               thermal plant types for new construction
                 CCStech(thermaltech)
         nucleartech(thermaltech)
         CCtech(thermaltech)
         renewtech(tech)                 renewable plant types for new construction
                 windtech(renewtech)
                 offshorewindtech(windtech)
                 solartech(renewtech)
*Storage
         storagetech(tech)               storage plant types for new construction
                 nonstoragetech(tech)                    non storage techs
                 ststoragetech(storagetech)
                 ltstoragetech(storagetech)
*CO2 removal
                 dacstech(tech)
                 notdacstech(tech)
*CANDIDATE TRANSMISSION LINES FOR CONSTRUCTION
*               ltech
*HOURS
         prmH(h)                        hour with peak net demand
         nonInitH(h)                   hours that are not initial horus of a block
*Building upgrade option
        bldgUp
*Interites
        interties
         ;

Alias (h,hh);
Alias (l,ll);

Parameters
*MAX NEW UNITS TO BE BUILT
        pNMaxRE(renewtech)
        pNMaxNuclear
        pNMaxCC
        pNMaxCCS(CCStech)
        pNMaxDACS(dacstech)
        pPMaxSto(storagetech)
        pEMaxSto(storagetech)
        pPERatio(ststoragetech)
*        pNMaxRegionalWind(z)
*        pNMaxRegionalSolar(z)
*        pNMaxTotalWind
*        pNMaxTotalSolar
*MINIMUM WIND AND SOLAR GENERATION
        pWSGenAsFracOfDemand
        pWindGenAsFracOfDemand
        pOSWCapReq
        pOSWExistingCap
*CAPACITY AND CAPACITY DERATING PARAMETERS 
         pCapactech(tech)                nameplate capacity of new builds for cost calculations [GW]
         pCapacderatetech(tech,h)       fractional derating of capacity 
*HEAT RATES [MMBtu/GWh]
         pHrtech(tech)                   heat rate of new builds [MMBtu per GWh]
*COST PARAMETERS
         pOpcosttech(tech)               total operational cost [thousandUSD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
         pFom(tech)                      fixed O&M cost [thousand$ per GW per yr]
         pOcc(tech)                      overnight capital cost [thousandUSD per GW]
                 pPowOcc(storagetech)                   occ for power capcity for lt storage
                 pEneOcc(storagetech)                   occ for energy capacity for lt storage
*                pCnse                                                  cost of nonserved energy [thousandUSD per GW]
*RAMP RATES [GW/hr]
         pRampratetech(tech)             up and down ramp rate of EGU assumed to be the same up & down [GW per hr]
*STORAGE PARAMETERS
         pEfficiencytech(storagetech)        round trip storage efficiency
         pChargeDischargeCapRatio       "ratio of charging to discharging" /1/
*EMISSIONS RATES [short ton/MMBtu]
         pCO2emratetech(tech)            CO2 emissions rate of potential new generators [short ton per MMBtu]
*EMISSIONS CAP AND COST
         pCO2emcap                       CO2 annual emissions cap [short tons]
*RENEWABLE GENERATION CAPS
         pCf(renewtech,h)                hourly capacity factors for potential new renewables
*FINANCIAL PARAMETERS
         pR                              discount rate
         pLife(tech)                     lifetime of tech [years]
                 pLifeline
         pCrf(tech)                      capital recovery factor
                 pCrfline
*ZONAL PARAMETERS
                pGenzonetech(tech)
                pLinecost(l)
         pNMaxLine(l)
*HOURLY RESERVE REQUIREMENTS [GW]
         pRegUpReqIncRE(renewtech,h)
                 pFlexReqIncRE(renewtech,h)
*RESERVE PROVISION PARAMETERS
         pRegeligibletech(tech)          candidate plant types eligible to provide reg reserves [0] or not [1]
         pFlexeligibletech(tech)
         pConteligibletech(tech)
         pMaxflexoffertech(tech)
         pMaxcontoffertech(tech)
         pMaxregupoffertech(tech)
*PLANNING RESERVE AND CONTRIBUTION PARAMETERS
         pPlanningreserve                planning margin reserve capacity [GW]
         pFOR(egu,h)                forced outage rate in each hour of existing gens
         pFORwind
         pFORsolar
         pFORtech(tech,h)           forced outage rate in each hour of candidate techs
         pPRMeligwindsolar          binary param indicating whether existing wind and solar can (1) or cannot (0) contribute to PRM
*BUILDING UPGRADE PARAMETERS
        pBldgUpgradeSavings(bldgUp,z,h)
        pBldgUpgradeCost(bldgUp)
*INTERTIE PARAMETERS
        pIntertieZone(interties)
        pIntertieCapacity(interties)
        pIntertiePrice(interties)
         ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load ststorageegu, ltstorageegu, tech, dacstech, thermaltech, CCStech, nucleartech, CCtech, renewtech, solartech, windtech, offshorewindtech, storagetech, ststoragetech, ltstoragetech, prmH
$load pNMaxRE, pNMaxNuclear, pNMaxCC, pNMaxCCS, pNMaxDACS, pPMaxSto, pEMaxSto
*$load pNMaxRegionalWind, pNMaxRegionalSolar, pNMaxTotalWind, pNMaxTotalSolar
$load pWSGenAsFracOfDemand, pWindGenAsFracOfDemand, pOSWCapReq, pOSWExistingCap
$load pCapactech, pCapacderatetech, pHrtech, pOpcosttech
$load pFom, pOcc, pPowOcc, pEneOcc, pRampratetech, pCO2emratetech
$load pEfficiencytech
$load pCO2emcap, pCf, pR, pLife, pLifeline
$load pGenzonetech, pLinecost, pNMaxLine
$load pRegUpReqIncRE, pFlexReqIncRE
$load pFlexeligibletech, pConteligibletech, pRegeligibletech
$load pPlanningreserve, pFOR, pFORwind, pFORsolar, pFORtech, pPRMeligwindsolar
$load bldgUp,pBldgUpgradeSavings,pBldgUpgradeCost
$load interties,pIntertieZone,pIntertieCapacity,pIntertiePrice
$gdxin

*DEFINE SET EXCLUSIONS
nonrenewegu(egu) = not renewegu(egu);
nonstoragetech(tech) = not storagetech(tech);
notdacstech(tech) = not dacstech(tech);
gentechs(tech) = not (dacstech(tech) + storagetech(tech));
*CALCULATE CAPITAL RECOVERY FACTOR
pCrf(tech) = pR / (1 - (1 / ( (1 + pR)**pLife(tech))));
pCrfline = pR / (1 - (1 / ( (1 + pR)**pLifeline)));
*CALCULATE PE RATIO FOR SHORTTERM STORAGE
pPERatio(ststoragetech) = pEMaxSto(ststoragetech)/pPMaxSto(ststoragetech);
*CALCULATE MAX RESERVE OFFERS
pMaxflexoffertech(tech) = pFlexeligibletech(tech)*pRampratetech(tech)*pRampratetoflexreservescalar;
pMaxcontoffertech(tech) = pConteligibletech(tech)*pRampratetech(tech)*pRampratetocontreservescalar;
pMaxregupoffertech(tech) = pRegeligibletech(tech)*pRampratetech(tech)*pRampratetoregreservescalar;

Variable
                 vVarcosttech(tech,h)
         vVarcostannual                             total variable costs for new and existing plants = variable O&M + fuel + emission costs [thousandUSD per yr]
                 vFixedcostannual                             total investment costs for new plants = fixed O&M + overnight capital costs [thousandUSD per yr]
         vGentech(tech,h)                hourly electricity generation by new plants [GWh]
*Emission variables
                 vCO2emstech(tech,h)
         vCO2emsannual                   co2 emissions in entire year from new and existing plants [short ton]
         ;

Positive variables
*Reserve requirements
         vRegupreserve(z,h)                amount of reg up reserves [GW]
         vFlexreserve(z,h)
*Generation and reserve variables
         vReguptech(tech,h)              hourly reg up reserves provided by new plants [GWh]
         vFlextech(tech,h)
         vConttech(tech,h)
*Storage
                 vStateofchargetech(storagetech,h)            "energy stored in storage unit at end of hour h (GWh)"
                 vChargetech(storagetech,h)                           "charged energy by storage unit in hour h (GWh)"
                 vPowBuiltSto(storagetech)                      built power capacity for storage
                 vEneBuiltSto(storagetech)                      built energy capacity for storage
*Line builds and flows
                vNl(l)
                vLinecapacnew(l)
*               vLinenewflow(ltech,h)
*Intertie flows
        vImport(interties,h)
                   ;
                   
Binary variables
        vNBldgUpgrade(bldgUp)
        ;

$include CEBuildVariable.gms
$include CEBldgUpgradeFix.gms

Equations
*Costs
                calcvarcoststech(tech,h)
                investmentcost
*Meet demand and reserves
         meetdemand(z,h)                   meet supply with demand
         meetreservemargin               meet planning reserve requirement with installed capacity
         setflexreserve(z,h)               determine quantity of required spin reserves
         setregupreserve(z,h)              determine quantity of required reg up reserves
         meetflexreserve(z,h)              meet spin reserve requirement
         meetcontreserve(z,h)              meet contingency reserve requirement
         meetregupreserve(z,h)             meet reg up reserve requirement
         limitshiftingdemandupper(z,h)          limit the amount of demand that can be shifted
         limitshiftingdemandlower(z,h)          limit the amount of demand that can be shifted
         meetshiftingdemand(z,h)
                limitLineFlows(l,h)
                linecapacnew(l,ll)

*Maximum build constraints
                 maxSolar(solartech)
                 maxWind(windtech)
                 maxNuclear(nucleartech)
                 maxCCS(CCStech)
                 maxCC(CCtech)
                 maxDAC(dacstech)
                 maxPSto(storagetech)
                 setVNSto(storagetech)
                 maxELTSto(ltstoragetech)
                 maxESTSto(ststoragetech)
                 maxL(l)
*                 maxTotalWind
*                 maxTotalSolar
*                 maxRegionalWind(z)
*                 maxRegionalSolar(z)
*Minimum wind & solar build constraint
                 lowerBoundWindSolarGenAsFracOfDemand
                 lowerBoundWindGenAsFracOfDemand
                 lowerBoundOffshoreWindCapacity
*Carbon emissions
                calcco2emstech(tech,h)
*Limit building upgrades to 1 upgrade
        limitBuildingUpgrades
*Interites
        limitImports(interties,h)
         ;

*****************CALCULATE COSTS*******************
*Variable costs for new techs
calcvarcoststech(tech,h).. vVarcosttech(tech,h) =e= vGentech(tech,h)*pOpcosttech(tech);

*Fixed costs = annual fixed O&M + fixed annualized capital costs
investmentcost..         vFixedcostannual =e= sum(nonstoragetech,vN(nonstoragetech)*pCapactech(nonstoragetech)*(pFom(nonstoragetech)+pOcc(nonstoragetech)*pCrf(nonstoragetech)))
                                                 + sum(storagetech,vPowBuiltSto(storagetech)*pPowOcc(storagetech)*pCrf(storagetech)+vEneBuiltSto(storagetech)*pEneOcc(storagetech)*pCrf(storagetech))
                                                                                                 + sum(l,vNl(l)*pLinecost(l)*pCrfline)
                                                                                                 + sum(bldgUp,vNBldgUpgrade(bldgUp)*pBldgUpgradeCost(bldgUp));
***************************************************

******************SYSTEM-WIDE GENERATION AND RESERVE CONSTRAINTS*******************
*Demand = generation by new and existing plants
meetdemand(z,h)..          sum(tech$[pGenzonetech(tech)=ORD(z)],vGentech(tech,h)) + sum(egu$[pGenzone(egu)=ORD(z)],vGen(egu,h)) +
                                                        sum(l$[pLinesink(l)=ORD(z)],vLineflow(l,h))*pTransEff + sum(interties$[pIntertieZone(interties)=ORD(z)],vImport(interties,h)) =g= (pDemand(z,h) + sum(bldgUp,vNBldgUpgrade(bldgUp)*pBldgUpgradeSavings(bldgUp,z,h)) + vShiftedDemand(z,h)
                                                        + sum(storageegu$[pGenzone(storageegu)=ORD(z)],vCharge(storageegu,h)) + sum(storagetech$[pGenzonetech(storagetech)=ORD(z)],vChargetech(storagetech,h))
                                                                                                                + sum(l$[pLinesource(l)=ORD(z)],vLineflow(l,h)));

*meetdemand(z,h)..          sum(tech$[pGenzonetech(tech)=ORD(z)],vGentech(tech,h)) + sum(egu$[pGenzone(egu)=ORD(z)],vGen(egu,h)) =g= pDemand(z,h);

*Demand response
limitshiftingdemandupper(z,h)..       vShiftedDemand(z,h) =l= pDemand(z,h)*pDemandShifter;
limitshiftingdemandlower(z,h)..       vShiftedDemand(z,h) =g= -1*pDemand(z,h)*pDemandShifter;
meetshiftingdemand(z,h)$(mod(ord(h),pDemandShiftingBlock)=1).. sum(hh$((ord(hh)>=ord(h) and (ord(hh)<=ord(h)+pDemandShiftingBlock-1))), vShiftedDemand(z,h))=e= 0;

*Meet planning reserve margin
meetreservemargin..       sum((thermaltech,prmH),pCapactech(thermaltech) * vN(thermaltech) * (1-pFORtech(thermaltech,prmH)) * (1-pCapacderatetech(thermaltech,prmH)))
                           + sum((storagetech,prmH),vPowBuiltSto(storagetech) * (1-pFORtech(storagetech,prmH)))
                           + sum((renewtech,prmH),pCapactech(renewtech) * vN(renewtech) * pCf(renewtech,prmH) * pFORtech(renewtech,prmH))
                           + sum((nonrenewegu,prmH),pCapac(nonrenewegu) * (1-pFOR(nonrenewegu,prmH)) * (1-pCapacderate(nonrenewegu,prmH)))
                           + sum((z,prmH),(pMaxgenwind(z,prmH) * pFORwind + pMaxgensolar(z,prmH) * pFORsolar) * pPRMeligwindsolar) =g= pPlanningreserve;
*                           + sum((bldgUp,z,prmH),vNBldgUpgrade(bldgUp)*pBldgUpgradeSavings(bldgUp,z,prmH));

*Define spinning and reg reserve requirement based on new builds
setflexreserve(z,h)..      sum(renewtech$[pGenzonetech(renewtech)=ORD(z)],vN(renewtech)*pCapactech(renewtech)*pFlexReqIncRE(renewtech,h)) + pFlexreserves(z,h) =e= vFlexreserve(z,h);
setregupreserve(z,h)..     sum(renewtech$[pGenzonetech(renewtech)=ORD(z)],vN(renewtech)*pCapactech(renewtech)*pRegUpReqIncRE(renewtech,h)) + pRegupreserves(z,h) =e= vRegupreserve(z,h);

*Meet spinning and regulation reserve requirements
meetflexreserve(z,h)..     sum(tech$[pGenzonetech(tech)=ORD(z)],vFlextech(tech,h)) + sum(egu$[pGenzone(egu)=ORD(z)],vFlex(egu,h)) =g= vFlexreserve(z,h);
meetcontreserve(z,h)..     sum(tech$[pGenzonetech(tech)=ORD(z)],vConttech(tech,h)) + sum(egu$[pGenzone(egu)=ORD(z)],vCont(egu,h)) =g= pContreserves(z,h);
meetregupreserve(z,h)..    sum(tech$[pGenzonetech(tech)=ORD(z)],vReguptech(tech,h)) + sum(egu$[pGenzone(egu)=ORD(z)],vRegup(egu,h)) =g= vRegupreserve(z,h);

*Limit line flows on new lines
linecapacnew(l,ll)$(pLinesource(l)=pLinesink(ll) and pLinesource(ll)=pLinesink(l)).. vLinecapacnew(l) =e= vNl(l) + vNl(ll);
limitLineFlows(l,h).. (pLinecapac(l)+vLinecapacnew(l))*(1-pLinecapacderate(l,h)) =g= vLineflow(l,h);

*Limit building upgrade investments to 1 set of packages
limitBuildingUpgrades .. sum(bldgUp,vNBldgUpgrade(bldgUp)) =l= 1;

*Limit imports on intertie
limitImports(interties,h) .. vImport(interties,h) =l= pIntertieCapacity(interties);
***********************************************************************************

*************UPPER AND LOWER GENERATION BOUNDS*****
*Set lower bound to zero for generation by storage and generating techs.
vGentech.lo(notdacstech,h) = 0;
*Set upper bound to zero for DACS techs. DACS max capacity is negative and has negative vGen values.
vGentech.up(dacstech,h) = 0;
***************************************************

******************BUILD DECISIONS******************
*Limit number builds to input value.
maxWind(windtech) .. vN(windtech) =l= pNMaxRE(windtech);
maxSolar(solartech) .. vN(solartech) =l= pNMaxRE(solartech);
maxNuclear(nucleartech) .. vN(nucleartech) =l= pNMaxNuclear;
maxCCS(CCStech) .. vN(CCStech) =l= pNMaxCCS(CCStech);
maxCC(CCtech) .. vN(CCtech) =l= pNMaxCC;
maxDAC(dacstech) .. vN(dacstech) =l= pNMaxDACS(dacstech);
maxPSto(storagetech) .. vPowBuiltSto(storagetech) =l= pPMaxSto(storagetech);
setVNSto(storagetech) .. vN(storagetech) =e= vPowBuiltSto(storagetech)/pCapactech(storagetech);
maxELTSto(ltstoragetech) .. vEneBuiltSto(ltstoragetech) =l= pEMaxSto(ltstoragetech);
maxESTSto(ststoragetech) .. vEneBuiltSto(ststoragetech) =e= pPERatio(ststoragetech)*vPowBuiltSto(ststoragetech);
maxL(l) .. vNl(l) =l= pNMaxLine(l);

*Constrain total and regional wind + solar deployment
*NEED TO MODIFY BELOW TO INCLUDE OFFSHORE WIND (STARTED BUT DIDN"T FINISH)
*maxTotalWind .. sum(windtech,vN(windtech)) =l= pNMaxTotalWind;
*maxTotalSolar .. sum(solartech,vN(solartech)) =l= pNMaxTotalSolar;
*maxRegionalWind(z).. sum(windtech$[pGenzonetech(windtech)=ORD(z)],vN(windtech)) =l= pNMaxRegionalOnshoreWind(z);
*maxRegionalWind(z).. sum(windtech$[pGenzonetech(windtech)=ORD(z)],vN(windtech)) =l= pNMaxRegionalOffshoreWind(z);
*maxRegionalSolar(z).. sum(solartech$[pGenzonetech(solartech)=ORD(z)],vN(solartech)) =l= pNMaxRegionalSolar(z);

*Set lower bound on wind and solar generation as percent of total demand
lowerBoundWindSolarGenAsFracOfDemand.. sum((renewegu,h),vGen(renewegu,h)) + sum((renewtech,h),vGentech(renewtech,h)) =g= sum((z,h),pWSGenAsFracOfDemand * pDemand(z,h));
lowerBoundWindGenAsFracOfDemand.. sum((windegu,h),vGen(windegu,h)) + sum((windtech,h),vGentech(windtech,h)) =g= sum((z,h),pWindGenAsFracOfDemand * pDemand(z,h));
lowerBoundOffshoreWindCapacity.. pOSWExistingCap + sum(offshorewindtech,vN(offshorewindtech)*pCapactech(offshorewindtech)) =g= pOSWCapReq;
***************************************************

******************RESERVE CONSTRAINTS******************
vFlextech.fx(tech,h)$[pMaxflexoffertech(tech)=0] = 0;
vConttech.fx(tech,h)$[pMaxcontoffertech(tech)=0] = 0;
vReguptech.fx(tech,h)$[pMaxregupoffertech(tech)=0] = 0;
*******************************************************

********CALCULATE CO2 EMISSIONS*************
calcco2emstech(tech,h) .. vCO2emstech(tech,h) =e= vGentech(tech,h)*pHrtech(tech)*pCO2emratetech(tech);
********************************************

