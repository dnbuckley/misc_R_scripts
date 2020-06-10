#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:42:12 2020

@author: dbuckley

First attempt to automate primer design by writing a script to interact
with zymo bisulfite primer seeker http://bps.zymoresearch.net/

Never written in python before so forgive the abysmal syntax
"""
import sys
import csv
import time
import re
import bs4 as bs
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import Select

# chromedriver should match current version of chrome installed
# may need to be updated from https://chromedriver.chromium.org/downloads
# if this is to be given to other people a better method will need to
# be devised...
def driver_init():
    driver = webdriver.Chrome(executable_path = "/Users/dbuckley/Desktop/bin/python_packages/chromedriver")
    driver.get("http://bps.zymoresearch.net/")
    driver.wait = WebDriverWait(driver, 30)
    return driver

def page_init(driver, seq,
              minAmpLen, maxAmpLen,
              minMeltTemp, maxMeltTemp,
              primerMinLen, primerMaxLen):
    driver.refresh()
    # set type based on seq.size
    tiling = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[4]/div[1]/div/input')
    regular = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[4]/div[2]/div/input')
    if (len(seq) <= 250):
        regular.click()
    else:
        tiling.click()
    # set defaults chris uses
    seqBox = driver.find_element_by_xpath('//*[@id="sequence"]')
    cpg_box = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[5]/div/input')
    cpg_box.click()
    submitButton = driver.find_element_by_xpath('//*[@id="submit"]')
    # min len
    minLenBox = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[3]/div[1]/input')
    minLenBox.clear()
    minLenBox.send_keys(minAmpLen)
    # max len
    maxLenBox = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[3]/div[2]/input')
    maxLenBox.clear()
    maxLenBox.send_keys(maxAmpLen)
    # min temp
    minTempBox = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[2]/div[1]/input')
    minTempBox.clear()
    minTempBox.send_keys(minMeltTemp)
    # max temp
    maxTempBox = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[2]/div[2]/input')
    maxTempBox.clear()
    maxTempBox.send_keys(maxMeltTemp)
    # primer min len
    primerMinLenBox = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[1]/div[1]/input')
    primerMinLenBox.clear()
    primerMinLenBox.send_keys(primerMinLen)
    # primer max len
    primerMaxLenBox = driver.find_element_by_xpath('//*[@id="form"]/div/div[2]/div[1]/div[2]/input')
    primerMaxLenBox.clear()
    primerMaxLenBox.send_keys(primerMaxLen)
    # enter sequence
    seqBox.clear()
    seqBox.send_keys(seq)
    # click the button
    submitButton.click()
    # wait for page to load
    driver.wait.until(EC.presence_of_element_located((By.XPATH, '//*[@id="table_length"]/label/select')))
    table_len = Select(driver.find_element_by_xpath('//*[@id="table_length"]/label/select'))
    table_len.select_by_visible_text('All')
    time.sleep(1)
    return

def parse_page(driver, out, ID, condition):
    html = driver.page_source
    # now that we have manipulated zymo we can extract the information
    # from the html
    soup = bs.BeautifulSoup(html, 'lxml') 
    table = soup.find('table', id = 'table')
    table_body = table.find('tbody')
    rows = table_body.find_all('tr')
    
    # the way they encoded this table is gnarly thus the weird parsing
    # and nested for loop
    failure = False
    for row in rows:
      sub_rows = row.find_all('td')
      amp = row.get('id')
      if amp is None:
        amp = "NA"
      forward = [ID, amp, strand, condition]
      reverse = [ID, amp, strand, condition]
      for sub in sub_rows:
        if sub.text == 'Show Amplicon':
          continue
        if sub.text == 'Sorry no valid primers found. Please try changing the parameters.':
          failure = True
          for i in range(0, 5):
            forward.append("NA")
            reverse.append("NA")
        else:
          f = sub.contents[0]
          r = sub.contents[2]
          forward.append(f)
          reverse.append(r)
      out.write('\t'.join(forward[:]) + '\n')    
      out.write('\t'.join(reverse[:]) + '\n')
    return failure
    

# print(str(sys.argv[1]))
seqs = open(str(sys.argv[1]))
# seqs = open("seqs.tsv")
n = str(sys.argv[2])
seqs = csv.reader(seqs, delimiter="\t")
out = open("zymo_output_" + n + ".tsv", "w+")
colNames = ["ID", "amp", "strand", "condition", "direction",
            "seq", "start", "size", "melting_temp"]
out.write('\t'.join(colNames[:]) + '\n')

# initialize zymo thing
driver = driver_init()

for r in seqs:
    seq = r[0]
    ID = r[1]
    strand = r[2]
    print("Processing", ID, "from", strand, "strand", "...")
    sys.stdout.flush()
    
    # Condition 1
    page_init(driver, seq,
              115, 180,
              60, 65,
              24, 36)
    fail1 = parse_page(driver, out, ID, "1")
    # Condition 2
    page_init(driver, seq,
              115, 180,
              55, 65,
              24, 36)
    fail2 = parse_page(driver, out, ID, "2")
    # Condition 3
    page_init(driver, seq,
              115, 225,
              55, 65,
              24, 36)
    fail3 = parse_page(driver, out, ID, "3")
    # Condition 4
    page_init(driver, seq,
              115, 225,
              55, 65,
              21, 36)
    fail4 = parse_page(driver, out, ID, "4")
    if (fail1 and fail2 and fail3 and fail4):
        print("No primers found under any condition.")
        sys.stdout.flush()
    
driver.close()
out.close()







