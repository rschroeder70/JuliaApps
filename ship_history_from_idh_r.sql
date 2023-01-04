--ship_history_from_idh_r.sql

-- WITH cteFRT AS (
--   SELECT xxtrf__site
--       ,xxtrf__ship_date
--       ,xxtrf__bol
--       ,xxtrf__mstr_bol
--      -- ,xxtrf__paid_date
--       --,xxtrf__scac_inv
--       --,xxtrf__scac_inv_date
--       ,SUM(xxtrf__chrg) xxtrf__chrg
--   FROM MFGPROPRD.dbo.xxtrf__chrg
--   GROUP BY xxtrf__site
--       ,xxtrf__ship_date
--       ,xxtrf__bol
--       ,xxtrf__mstr_bol
--      -- ,xxtrf__paid_date
--      -- ,xxtrf__scac_inv
--      -- ,xxtrf__scac_inv_date
-- )
SELECT YEAR(ih_ship_date)*100 + MONTH(ih_ship_date) ShipMonth, 
  CAST(ih_ship_date AS DATE) ShipDate, 
  idh_site ShipSite,
  --ih_ship,
  MS.MST,
  MS.MSTSetup,
  --cm_name ShipToName,
  --LEFT(cm_name, 5) ShipToName5,
  MS.MSTName,
  --ad_city City, 
  MS.MSTCity,
  --ad_state State,
  MS.MSTState, 
  --ad_zip Zip, 
  MS.MSTZip,
  --ad_country Country,
  MS.MSTCountry,
  ih_shipvia shipVia, 
  UPPER(ih_fr_terms) FrtTerms, 
  CASE WHEN UPPER(ih_fr_terms) = 'PREPAID' THEN 'PPD' ELSE 'CPU' END FrtTermsFix,
  --ih_nbr SalesOrder,
  --ih_bol BillOfLading, 
  LEFT(ih_bol, 9) MasterBOL, 
  --xxih__pro_num ProNbr, 
  xxih__custpu CustPUFlag, 
  --xxih__stop_nbr StopNbr,
  DENSE_RANK() over (partition by left(ih_bol, 9) order by MS.MST ASC) as StopNbrFix,
  -- DENSE_RANK() over (partition by left(ih_bol, 9) order by LEFT(cm_name, 5) + ad_zip) +
  --      DENSE_RANK() over (partition by left(ih_bol, 9) order by LEFT(cm_name, 5) + ad_zip DESC) - 1 as StopsFix,
  -- DENSE_RANK() over (partition by left(ih_bol, 9) order by ih_ship) +
  --      DENSE_RANK() over (partition by left(ih_bol, 9) order by ih_ship DESC) - 1 as BOLShipToCount,  
  DENSE_RANK() over (partition by left(ih_bol, 9) order by MS.MST) +
      DENSE_RANK() over (partition by left(ih_bol, 9) order by MS.MST DESC) - 1 as BOLMasterShipToCount,  
  xxih__load Load,
  xxih__truck Truck,
  CASE WHEN xxih__truck IN ('NLAT', 'QPMT', 'SCDS', 'SNCK', 'SNLU' , 'RBIN' , 'ANVI' , 'HUBG' , 'GRSF' , 'COYY' , 'BTIU', 'HJBI', 'HBGI') THEN 'IM'
            WHEN xxih__truck IN ('CNWY','CTII', 'CWCE', 'CWWE', 'DHRN', 'FXFE', 'HMES', 'PITD', 'RDWY', 'RETL', 'SAIA', 'UPFG', 'UPGF', 'YFSW' , 'XPO' , 
                                 'YFSY', 'FXNL' , 'EXLA' , 'SEFL' , 'PYLE' , 'DHRN' , 'LAXV' , 'TONI' , 'HJBA' , 'NPME') THEN 'LTL'
            WHEN xxih__truck IN ('CPU' , 'cpu' , 'CUST') THEN 'CPU'
            WHEN xxih__truck IN ('PAGE') THEN 'OCEAN' 
            WHEN xxih__truck IN ('FEAF','UPS', 'FEDX', 'DBYA', 'UPSA') THEN 'AIR' 
            WHEN xxih__truck IN ('AMNO', 'BIOS', 'CFAA', 'CRCR', 'GATI', 'HANY', 'HJBT', 'IWXP', 'KMTN', 'LDWY', 'LRGR', 'RBTW', 'RISI', 'RLOG', 
                                 'SCNN', 'SWFT', 'TBIF', 'USIT', 'WMSN', 'WSXI', 'SNCY', 'MKTN', 'CLLQ', 'WSGS', 'TRCR', 'GRTV', 'SDAT'  , 
                                 'USIT' , 'usit' , 'HJBB' , 'WENP' , 'ANIQ' , 'HJBL' , 'AXLL', 'CLNC', 'CVYI', 'LFRN', 'HYBL', 'ECHS', 'EGBL',
                                 'ABXN', 'DBOI', 'AVGW', 'XPON', 'DETR', 'DCLK', 'KLFM', 'CMMS', 'CEQT', 'PTWT', 'TTMS', 'HNRF', 'BKLI', 'EOSA',
                                 'CHXD', 'XPOL', 'LTCA', 'CONV', 'PLYD', 'QRKQ', 'ERWN', 'SWIF', 'FTMG', 'CVYE', 'VITI', 'XPOX', 'NAAF', 'WOLH',
                                 'RCXV', 'TTFO', 'NWLS', 'GULF', 'PTWB', 'HGTN', 'DDOG', 'TFFO', 'GRLO', 'DEBE', 'CNVY', 'JBHU', 'PAPT') THEN 'TL'
            ELSE 'UNKNOWN' END ModeFix,

  pt_part_type PartType,
  pt_group PrintGroup,
  SUM(ISNULL(idh_qty_ship, 0)) Cases,
  SUM(ISNULL(idh_qty_ship, 0) * xxpt__case_pk * 0.001) MUnits,
  SUM(ISNULL(idh_qty_ship, 0) * pt_size) Cube,
  SUM(ISNULL(idh_qty_ship/P.CasepPallet, 0)) as Pallets,
  ISNULL(xxih__base_chrg, 0) xxih__base_chrg,
  ISNULL(FRT.Miles, 0) Miles,
  ISNULL(FRT.Freight, 0) Freight, ISNULL(FRT.EstimatedFreight, 0) EstimatedFreight, 
  ISNULL(FRT.InvFrt, 0) InvFrt, ISNULL(FRT.ActualFreight, 0)ActualFreight, ISNULL(FRT.PaidFlag, 0) PaidFlag
 -- xxih__inv_tot gross_sales
FROM MFGPROPRD.dbo.ih_hist_all 
  LEFT JOIN MFGPROPRD.dbo.idh_hist_all ON ih_domain = idh_domain AND ih_nbr = idh_nbr
  LEFT JOIN MFGPROPRD.dbo.pt_mstr_all ON idh_domain = pt_domain AND idh_part = pt_part 
  LEFT JOIN MFGPROPRD.dbo.ad_mstr_all ON ad_mstr_all.ad_addr = ih_ship AND ad_mstr_all.ad_domain = '037'
  LEFT JOIN MFGPROPRD.dbo.cm_mstr_all as cmst ON ad_mstr_all.ad_addr = cmst.cm_addr AND cmst.cm_domain = '037' 
    AND ad_mstr_all.ad_domain = cmst.cm_domain
--  LEFT JOIN MFGPROPRD.dbo.ptp_det ON ptp_det.ptp_site = idh_site 
--    AND ptp_part =  idh_part
  LEFT JOIN dbo.Pallets P on P.Part = idh_part AND P.Site = idh_site
  LEFT JOIN dbo.Freight FRT ON FRT.ShipSite = ih_site AND FRT.SalesOrderNo = ih_nbr AND FRT.ShipDate = ih_ship_date
  LEFT JOIN dbo.MasterShipTo MS ON MS.ST = ih_ship
  --ON m.xxtrf__domain = dbo.ih_hist_all.ih_domain AND dbo.ih_hist_all.ih_nbr = m.xxtrf__so_nbr AND dbo.ih_hist_all.ih_ship_date = m.xxtrf__ship_date 
--  LEFT JOIN cteFRT FRT ON xxtrf__site = ih_site AND xxtrf__bol = ih_bol AND xxtrf__ship_date = ih_ship_date
WHERE ih_domain = '037'
       --AND ih_channel  in ('IN','DM','RE','II','IC','DI')
  AND xxih__chr01 in ('INV')
  AND ih_ship_date >= dbo.FirstDayofCurrentYear()
  AND cm_name NOT LIKE '%SAMPL%'
  AND idh_site != '531'
  --AND Left(ih_bol, 9) in ('045-56996', '031-96291')  -- Test of window functions; 1st is 1 stop, second is 2
  --AND left(ih_bol, 9) = '031-30884'
GROUP BY ih_ship_date, 
  idh_site, 
  --ih_nbr,
  --ih_ship,
  MS.MST,
  --cm_name,
  --ad_city, ad_state, ad_zip, ad_country,
  MS.MSTCity, MS.MSTCountry, MS.MSTName, MS.MSTSetup, MS.MSTState, MS.MSTZip,
  ih_shipvia, ih_fr_terms,
  pt_group, 
  pt_part_type,
  --ih_bol,
  left(ih_bol, 9),
  xxih__load,
  xxih__truck, 
  --xxih__pro_num, 
  --xxih__stop_nbr, 
  xxih__custpu,
  xxih__base_chrg, 
  FRT.Miles,
  FRT.Freight, FRT.InvFrt, FRT.ActualFreight, FRT.EstimatedFreight, FRT.PaidFlag

/*

select idh_site, YEAR(ih_ship_date), MONTH(ih_ship_date), xxih__chr01, ih_fr_terms,
SUM(idh_qty_ship * xxpt__case_pk * 0.001) MUnits
FROM MFGPROPRD.dbo.ih_hist_all 
  LEFT JOIN MFGPROPRD.dbo.idh_hist_all ON ih_domain = idh_domain AND ih_nbr = idh_nbr
  LEFT JOIN MFGPROPRD.dbo.pt_mstr_all ON idh_domain = pt_domain AND idh_part = pt_part 
where idh_site = '045' 
and ih_domain = '037'
AND xxih__chr01 in ('INV')
AND YEAR(ih_ship_date) = '2019'
AND MONTH(ih_ship_date) = '1'
GROUP BY idh_site, YEAR(ih_ship_date), MONTH(ih_ship_date), xxih__chr01,
ih_fr_terms


*/