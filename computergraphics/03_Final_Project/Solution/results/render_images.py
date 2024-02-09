import subprocess
import json

n = 3

def render_image(samples : int, randomizer : int):
    time = .0
    for i in range(n):
        output = subprocess.run(
            f'/usr/bin/time -f "{{\\"wall\\": %e}}" ./bin/PathTracing {samples} {randomizer}',
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE
        ).stderr.decode()
        time += json.loads(output)["wall"]
    time /= n
    print(time)

    return time

def calculate_psnr(samples :int, randomizer :int):
    ratio = subprocess.run(
        f'./bin/Image_diff {samples} {randomizer}',
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE
    )

    return json.loads(ratio.stderr.decode())["psnr"]

times = {}
psnr = {}
randomizers = ["drand48", "halton23", "stratified", "blue noise"]

for idx, rand in enumerate(randomizers):
    times[rand] = {}
    psnr[rand] = {}
    for samples in range(1, 116, 5):
        # print(f'[*] Rendering image with {samples} samples with random sampler {rand}...')
        # times[rand][samples] = render_image(samples, idx)
        # print("[+] Done")
        psnr[rand][samples] = calculate_psnr(samples, idx)

# with open('data/render.times', 'w+') as outfile:
#     outfile.write(json.dumps(times, indent=4))

with open('data/render.psnr', 'w+') as outfile:
    outfile.write(json.dumps(psnr, indent=4))
