from flask import Flask, request, jsonify
from BubbleFlask import BubblePoint

app = Flask(__name__)
app.config["DEBUG"] = True


@app.route('/', methods=['GET','POST'])
def home():
    if request.method == 'POST':
        data = request.json
        bubblepoint = BubblePoint(data)
        return bubblepoint

    elif request.method == 'GET':
        return "<h1>Distant Reading Archive</h1><p>This site is a prototype API for distant reading of science fiction novels.</p>"

app.run()