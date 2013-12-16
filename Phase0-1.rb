require 'rubygems'
require 'bundler/setup'

require 'cytogenetics'
require 'oracle_query'
require 'logger'
require 'yaml'

# Cytogenetics Gem log
log = Logger.new(STDOUT)
log.level = Logger::INFO
Cytogenetics.logger = log

# Connect to local Oracle DB containing a copy of the Mitelman DB
query = Oracle.new("system", "manager", "localhost:1521/XE")

#
# Phase 0 - Filtering
#

# Fetch list of karyotype samples
#
# -> Account for KARYSHORT/KARYLONG columns (karyotypes longer than 255bytes are added
#    to KARYLONG while KARYSHORT retains the first 255bytes)
#
# -> Removal of all karyotypes containing "idem" and "/" (cannot currently be parsed
#    by CytogeneticsParser)
#
# -> Removal of all karyotypes containing unknown/unclear aberrations ("?")
#
# TIP : Add 'WHERE ROWNUM <= 5000' before 'ORDER BY' to limit size of result set

query.select("SELECT * FROM (
(SELECT REFNO, CASENO, KARYLONG FROM SYSTEM.CYTOGENINV
   WHERE KARYLONG IS NOT NULL AND NOT (KARYLONG LIKE '%idem%' OR KARYLONG LIKE '%/%' OR KARYLONG LIKE '%?%'))
UNION ALL
(SELECT REFNO, CASENO, KARYSHORT FROM SYSTEM.CYTOGENINV
   WHERE KARYLONG IS NULL AND NOT (KARYSHORT LIKE '%idem%' OR KARYSHORT LIKE '%/%' OR KARYSHORT LIKE '%?%')))
WHERE ROWNUM <= 60000
ORDER BY REFNO,CASENO", "karylist_P00")

query.drop("table SYSTEM.CYTOGENINV_processed", "table")
query.create("SYSTEM.CYTOGENINV_processed", [["REFNO", 'NUMBER'],["CASENO", 'varchar2 (14 BYTE)'],
                                             ["CHR_GAIN",'NUMBER'],["CHR_LOSS",'NUMBER'],["TRANS",'NUMBER'],["TAIL_DEL",'NUMBER'],
                                             ["TYPE",'varchar2 (200 BYTE)'],["MORPH",'NUMBER']])

# Array of parsable karyotypes (no errors when fed into CytogenicsParser)
karylist_P01 = []

# Hash containing Number of rearrangement events
events_profile = {:CHR_GAIN => 0,
                  :CHR_LOSS => 0,
                  :TRANS => 0,
                  :TAIL_DEL => 0}

# Insert additional karyotypes for analysis here
#$karylist_P00.push("FOO")
#$karylist_P00 = [["0","0","47,XY,der(1)t(1;18)(p36;q21),t(14;18)(q32;q21),+der(18)t(12;18)(p11;q21),+der(18)t(14;18)(q32;q21)"]]

# Start karyotype analysis
$karylist_P00.each do |(refno,caseno,kary)|

    begin
        kt = Cytogenetics.karyotype(kary)

    rescue => error
        puts "Error parsing karyotype : #{error.class}"

    else

        events_profile[:CHR_GAIN] = kt.aberrations.values_at(:gain).flatten.compact.count
        events_profile[:CHR_LOSS] = kt.aberrations.values_at(:loss).flatten.compact.count
        events_profile[:TRANS] = kt.aberrations.values_at(:trans).flatten.compact.count
        events_profile[:TAIL_DEL] = kt.aberrations.values_at(:del).flatten.compact.count

        query.select("SELECT MORPH FROM SYSTEM.CYTOGEN WHERE CASENO = '#{caseno.to_i}' AND REFNO = #{refno.to_i}","morph")
        $morph = $morph.flatten.compact.first.to_i
        query.select("SELECT BENAMNING FROM SYSTEM.KODER WHERE KOD =  '#{$morph}'", "type")

        query.insert("SYSTEM.CYTOGENINV_processed", ["'#{refno.to_i}', '#{caseno.to_i}',
                                                      '#{events_profile[:CHR_GAIN]}',
                                                      '#{events_profile[:CHR_LOSS]}',
                                                      '#{events_profile[:TRANS]}',
                                                      '#{events_profile[:TAIL_DEL]}',
                                                      '#{$type.to_s}',
                                                      '#{$morph}'"], false )

        karylist_P01.push([refno,caseno,kary])

    end

end

# Phase completion information
puts "\nPhase 1 - Invalid karyotypes : #{$karylist_P00.count - karylist_P01.count}"
puts "Phase 1 - Final count : #{karylist_P01.count}"
