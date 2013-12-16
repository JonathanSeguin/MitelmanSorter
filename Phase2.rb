require 'rubygems'
require 'bundler/setup'

require 'bigdecimal'
require 'oracle_query'
require 'logger'
require 'yaml'

# Connect to local Oracle DB containing a copy of the Mitelman DB
query = Oracle.new("system", "manager", "localhost:1521/XE")

# Factorial function for integers
class Integer
  def fact
    (1..self).reduce(:*) || 1
  end
end

# Combination function based on the factorial function
def combi (a,b)
    if b > a then
        return BigDecimal('0.0')
    else
        return BigDecimal(BigDecimal(a.fact) / (BigDecimal(b.fact) * BigDecimal((a - b).fact)))
    end
end

# Hypergeometric score
def hypergeometric (n,m,n1,k)
  score = BigDecimal('0.0')
  # Summation
  for i in k..[n1,m].min do
    if combi(n,n1) != BigDecimal('0.0') then
        score = score + ((combi(m,i)*combi(n-m,n1-i)) / combi(n,n1))
    end
  end
  return score
end


#
# Phase 2 - Clustering
#

significant_pairs = 0

cancer_classes = [1103, 1107, 1108, 1203, 1501, 1505, 1509, 1702, 1810, 1820, 3003, 8804]
cancer_classes_pairs = cancer_classes.combination(2).to_a

cancer_classes_pairs.each do |(class_A, class_B)|
    query.select("SELECT *
                 FROM SYSTEM.CYTOGENINV_processed
                 WHERE MORPH = '#{class_A}'","c1")

    query.select("SELECT *
                 FROM SYSTEM.CYTOGENINV_processed
                 WHERE MORPH = '#{class_B}'","c2")

    omega = $c1 + $c2
    significant = 0

    puts "\n\n\n=== NEXT PAIRS ==="
    puts "Class A : #{$c1.first[6]}"
    puts "Class B : #{$c2.first[6]}"

    # Thresholds
    for t in 0..4 do

        d1_CHR_GAIN, d1_CHR_LOSS, d1_TRANS, d1_TAIL_DEL = [],[],[],[]

        omega.each do |kt|

            kt[2] <= t ? d1_CHR_GAIN.push(kt) : nil
            kt[3] <= t ? d1_CHR_LOSS.push(kt) : nil
            kt[4] <= t ? d1_TRANS.push(kt) : nil
            kt[5] <= t ? d1_TAIL_DEL.push(kt) : nil

        end

        puts "\n\nTHRESHOLD : #{t}"
        puts "\nCHR_GAIN"

        n = omega.count
        n1 = $c1.count
        m = d1_CHR_GAIN.count
        k = ($c1 & d1_CHR_GAIN).count

        p_value = BigDecimal('2.0')*[hypergeometric(n, n1, m, k), hypergeometric(n, n1, n - n1, m - k)].min

        if p_value < BigDecimal('0.001') && p_value != 0.0 then
            significant = 1
            puts "p-value = #{p_value}"
        end


        puts "\nCHR_LOSS"
        m = d1_CHR_LOSS.count
        k = ($c1 & d1_CHR_GAIN).count
        p_value = BigDecimal('2.0')*[hypergeometric(n, n1, m, k), hypergeometric(n, n1, n - n1, m - k)].min

        if p_value < BigDecimal('0.001') && p_value != 0.0 then
            significant = 1
            puts "p-value = #{p_value}"
        end


        puts "\nTRANS"
        m = d1_TRANS.count
        k = ($c1 & d1_CHR_GAIN).count
        p_value = BigDecimal('2.0')*[hypergeometric(n, n1, m, k), hypergeometric(n, n1, n - n1, m - k)].min

        if p_value < BigDecimal('0.001') && p_value != 0.0 then
            significant = 1
            puts "p-value = #{p_value}"
        end


        puts "\nTAIL_DEL"
        m = d1_TAIL_DEL.count
        k = ($c1 & d1_CHR_GAIN).count
        p_value = BigDecimal('2.0')*[hypergeometric(n, n1, m, k), hypergeometric(n, n1, n - n1, m - k)].min

        if p_value < BigDecimal('0.001') && p_value != 0.0 then
            significant = 1
            puts "p-value = #{p_value}"
        end

    end

    significant_pairs = significant_pairs + significant
end

puts "#{significant_pairs} / 66 significant"

