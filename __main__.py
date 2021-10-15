# ------------------------------------------------------------------------------
# Description:  Startet die Anwendung, andere Funktionen werden in der GUI auf-
#               gerufen 
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-04      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:      https://realpython.com/python-main-function/
# ------------------------------------------------------------------------------
# Imports:      
from gui import app     
# ------------------------------------------------------------------------------

def main():                                         # main() should contain code that you want to run when the file is executed
    app.start()

if __name__ == "__main__":                          # Calling main()-function when script is excecuted     
    main()                                         

# ------------------------------------------------------------------------------
# Aufbau des Codes:
# 1. main.py startet die Applikaion indem app.py aufgerufen wird
# 2. app.py startet die GUI und erstellt die Objekte aus building.py
# 3. Durch klicken eines Knopfes in der GUI ruft app.py die dem gewählten Trag-
#    werk entsprechende str_....py auf
# 4. In der str_....py werden die Tragwerksspezifischen Eigenschaften berechnet
#    und die Elemente mithilfe der Nachweise und Tragwerksübergreifender Berech-
#    nungen in calculations.py und fea.py bemessen und in den zuvor erstellten
#    Objekten gespeichert 