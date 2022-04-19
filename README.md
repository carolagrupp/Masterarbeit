# premiumForHeight
Programmcode zur parametrischen Analyse des Materialbedarfs von Aussteifungssystemen bei Hochhäusern



# structure

Aufbau des Codes:
1. main.py startet die Applikaion indem app.py aufgerufen wird
2. app.py startet die GUI und erstellt die Objekte aus building.py
3. Durch klicken eines Knopfes in der GUI ruft app.py die dem gewählten Tragwerk entsprechende str_....py auf
4. In der str_....py werden die Tragwerksspezifischen Eigenschaften berechnet und die Elemente mithilfe der Nachweise und Tragwerksübergreifender       Berechnungen in calculations.py und fea.py bemessen und in den zuvor erstellten Objekten gespeichert

In Wind Data ist die Vorbereitung zur Berücksichtigung von dynamischen Windlasten enthalten. Es können Windtunnelergebnisse ausgelesen und die Umwandlung der Daten findet statt. Noch ist es nicht möglich diese an dem System aufzubringen.


## examples
Es gibt die Möglichkeit Kerntragwerke, Rahmentragwerke, Framed Tube, Bundeld Tube und Outrigger zu analysieren.

## gui
Über die GUI können Gebäudegeometrie, Lasten und Parameter vorgegeben und variiert werden. 

# feastruct
Zur Berechnung des Systems als zweidimensionales finite Elemente System wurde das Skript feastruct hinterlegt und leicht erweitert. 
Zu Informationen zu den Rechten siehe feastruct/License



