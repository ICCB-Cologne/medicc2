import argparse
import os
import pickle


parser = argparse.ArgumentParser()
parser.add_argument("task_dir", 
                    help = "a path to the task directory")

args = parser.parse_args()

with open(os.path.join(args.task_dir, f'tasks.pickle'), 'rb') as f:
    task_idxs = pickle.load(f)

print(f"Found {len(task_idxs)} tasks.")

for idx in task_idxs:
    print(f"Running task {idx}")

    with open(os.path.join(args.task_dir, f'task_{idx}.pickle'), 'rb') as f:
        task = pickle.load(f)

    result = task[0](*task[1], **task[2])

    with open(os.path.join(args.task_dir, f'task_{idx}_result.pickle'), 'wb') as f:
        pickle.dump(result, f)

print(f"Finished tasks")
