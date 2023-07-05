from flask import Flask, request, render_template, send_file
from HCR_prober_generator import main


app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        seq = request.form['seq']
        probe_size = int(request.form['probe_size'])
        initiator_type = str(request.form['initiator_type'])
        polyN = int(request.form['polyN'])
        min_gc = float(request.form['min_gc'])
        max_gc = float(request.form['max_gc'])
        #blastdb = request.form['blastdb']
        output = 'output.csv'
        main(seq, probe_size, initiator_type, polyN, min_gc, max_gc, output)
        return send_file(output, as_attachment=True)
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
