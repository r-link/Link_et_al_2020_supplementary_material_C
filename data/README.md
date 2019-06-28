Meta information for the dataset
================

## Information

The present dataset is based on day-level aggregates of radial profiles
of *sap flux per section* (SFS) measured with the heat field deformation
(HFD) method on 38 trees belonging to 8 species of Costa Rican tropical
dry forest species measured in four campaigns of 5–7 days in October and
November 2015 in the Estación Experimental Forestal Horizontes,
Guanacaste, Costa Rica.

Detailed information about the sap flow meausrements measurement can be
found in the main text of the accompanying publication.

| Variable | Unit         | Example                 | Explanation                                                                                                                                                           |
| -------- | ------------ | ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| species  | \-           | *Tabebuia impetiginosa* | Tree species                                                                                                                                                          |
| tree     | \-           | H\_BS2\_01\_03          | Tree ID (Site\_Plot\_Subplot\_StemNo)                                                                                                                                 |
| campaign | \-           | 1                       | Measurement campaign (numbered in increasing order) - only 10 trees could be measured at a time due to the number of available sensors and batteries                  |
| date     | yyyy-mm-dd   | 2015-10-08              | Measurement day (data are aggregated to daily averages to average out diurnal changes)                                                                                |
| xylemrad | cm           | 10.4159                 | Approximate xylem radius at the height of sensor installation (calculated from circumference assuming a perfectly circular stem, and substracting the bark thickness) |
| depth    | cm           | 0.5                     | Sensor depth as distance from the cambium                                                                                                                             |
| reldepth | \-           | 0.1441                  | Sensor depth expressed as a fraction of the xylem radius                                                                                                              |
| flux     | cm³ cm⁻² s⁻¹ | 5.7868                  | Raw sap flux per section measured by the HFD method                                                                                                                   |
| WD       | g cm⁻³       | 0.6138                  | Sapwood wood density measured from samples excised with an increment corer (dry mass at 105°C over wet volume)                                                        |
| height   | m            | 12.65                   | Tree height determined from 4 consecutive measurements with an Haglöfs Vertex IV + Transponder                                                                        |
| ASI      | mm yr⁻¹      | 8.7306                  | Absolute annual stem diameter increment (slope of dendrometer readings from 2015-2017 against time)                                                                   |
