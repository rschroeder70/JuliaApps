module LoadForecast
export fcst_df

using ODBC
using DataFrames

conn = ODBC.Connection("CUPS", "crystal", "Bobject\$")

query = """ 
    SELECT ForecastSet,
        MfgSite,
        SiteName MfgSiteName,
        PrintGroup,
        PartType,
        SUM(OLMUnits) OLMUnits
    FROM dbo.FinanceForecastV1
        LEFT JOIN dbo.Sites ON Sites.SiteCode = MfgSite
    WHERE ForecastSet = '202205'
        AND YearMonth BETWEEN '202206' AND '202305'
        --AND PrintGroup = 'COKECORP'
    GROUP BY ForecastSet,
        MfgSite, 
        SiteName,
        PrintGroup,
        PartType;
"""
fcst_df = DBInterface.execute(conn, query) |> DataFrames.DataFrame

DBInterface.close!(conn)

end