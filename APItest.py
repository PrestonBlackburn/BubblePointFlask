import requests

url = 'http://127.0.0.1:5000/'
headers = {'Content-Type':'Application/json'}
data = '{"Temperature":74,"Hydrogen Sulfide":0,"Nitrogen":0.019264,"Carbon Dioxide":0.176679,"Methane":2.216492,"Ethane":2.131195,"Propane":2.672874,"Isobutane":1.531791,"n-Butane":3.01741,"n-Pentane":5.660723,"n-Hexane":8.777826,"n-Heptane":73.795745,"n-Octane":0,"n-Nonane":0,"n-Decane":0,"Benzene":0,"Toluene":0,"Xylenes":0,"Cyclohexane":0}'

r = requests.post(url, headers = headers, data = data)

print(r)
print(r.text)