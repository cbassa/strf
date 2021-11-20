#!/usr/bin/bash

DB_API_TOKEN="4f20a493a3f5fd85074d61db944365bf94987613"
URL="https://db-satnogs.freetls.fastly.net/media/artifacts/b4975058-04eb-4ab7-9c40-9bcce76d94db.h5"
echo "Authorization: Token $DB_API_TOKEN"
curl -H "Authorization: Token $DB_API_TOKEN" "$URL"

