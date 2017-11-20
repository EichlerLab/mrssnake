import os
import boto3
from subprocess import check_output

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dynamo_table", help="Name of dynamodb table to create")
    parser.add_argument("key_name", help="Name of dynamodb key")
    parser.add_argument("input_table", help="Newline-delimited list of sample names")

    parser.add_argument("--region", default="us-east-1", help="Region to create table in (Default: %(default)s)")
    args = parser.parse_args()

    client = boto3.client('dynamodb',
                          region_name=args.region,
                          aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
                          aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"]
                         )
    dynamodb = boto3.resource('dynamodb',
                              region_name=args.region,
                              aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
                              aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"]
                             )

    table = dynamodb.Table(args.dynamo_table)
    try:
        table.item_count
        print("Table already exists. Exitting...")
        sys.exit(1)
    except client.exceptions.ResourceNotFoundException:
        table = dynamodb.create_table(
            TableName=args.dynamo_table,
            AttributeDefinitions=[
                {
                    'AttributeName': 'SimonsID',
                    'AttributeType': 'S'
                },
            ],

            KeySchema=[
                {
                    'AttributeName': 'SimonsID',
                    'KeyType': 'HASH'
                },
            ],
            ProvisionedThroughput={
                'ReadCapacityUnits': 5,
                'WriteCapacityUnits': 5
            }
        )
    table.meta.client.get_waiter('table_exists').wait(TableName=args.dynamo_table)

    samples = []
    with open(args.input_table, "r") as input_samples:
        for line in input_samples:
            samples.append(line.rstrip())

    samples = list(set(samples))
    for sample in samples:
        table.put_item(
            Item={
                args.key_name: sample,
            }
        )
