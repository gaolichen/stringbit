@echo off
for %%f in (*.tex) do pdflatex %%f
del *.aux
del *.log
@echo on