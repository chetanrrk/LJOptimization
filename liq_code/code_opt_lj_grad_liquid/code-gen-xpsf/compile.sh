#!/bin/bash

g++ -O2 -o gen_xpsf gen_xpsf.cpp -static

cp gen_xpsf ..
chmod g+rx ../gen_xpsf
