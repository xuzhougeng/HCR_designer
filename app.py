from flask import Flask, request, render_template, send_file
from HCR_prober_generator import main


app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        name = request.form['name']
        seq = request.form['seq']
        probe_size = int(request.form['probe_size'])
        initiator_type = str(request.form['initiator'])
        polyN = int(request.form['polyN'])
        min_gc = float(request.form['min_gc'])
        max_gc = float(request.form['max_gc'])
        
        
        blast = request.form['blast'] == "y"
        blastdb = None
        if blast:
            blastdb = request.form['blastdb']
            if blastdb == 'Arabidopsis':
                balstdb = 'db/Arabidopsis'

        output = 'output.csv'
        
        main(name,seq, probe_size, initiator_type, polyN, min_gc, max_gc, output, blastdb)
        return send_file(output, as_attachment=True)
    
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
