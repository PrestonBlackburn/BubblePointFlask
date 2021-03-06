from flask import Flask, request, jsonify
from BubbleFlask import BubblePoint

app = Flask(__name__)

@app.route('/')
def home():
    return "<h1>Distant Reading Archive</h1><p>This site is a prototype API.</p>"


@app.route('/API', methods=['POST'])
def API():

    data = request.json
    bubblepoint = BubblePoint(data)
    return bubblepoint
    
if __name__ == "__main__":
    app.run()
